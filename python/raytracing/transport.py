"""Coherent branching transport over explicit material regions."""
from dataclasses import dataclass, field
import math
import numpy as np
from .materials import fresnel
from .scene import curvature, edge_distance, unit
from .beams import State, Segment, frame, propagate


@dataclass
class SolverConfig:
    frequency_hz: float = 77e9
    max_interactions: int = 16
    max_distance: float = 20.0
    min_power_fraction: float = 1e-10
    max_branches: int = 100000
    max_loss_ratio: float = 0.01

    def validate(self):
        if not math.isfinite(self.frequency_hz) or self.frequency_hz <= 0 or not math.isfinite(self.max_distance) or self.max_distance <= 0:
            raise ValueError("Frequency and max distance must be positive")
        if not isinstance(self.max_interactions,int) or self.max_interactions<0 or not isinstance(self.max_branches,int) or self.max_branches<1:
            raise ValueError("Invalid interaction/branch budget")
        if not 0 <= self.min_power_fraction < 1 or not 0 <= self.max_loss_ratio <= .01:
            raise ValueError("Require 0 <= power threshold < 1 and 0 <= weak-loss ratio <= 0.01")


@dataclass
class Terminal:
    branch: int
    reason: str
    position: np.ndarray
    power: float
    history: tuple


@dataclass
class Simulation:
    config: SolverConfig
    scene: object
    segments: list = field(default_factory=list)
    terminals: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    launched_power: float = 0.0
    absorbed_power: float = 0.0
    interface_flux_residual: float = 0.0

    def field_at(self, position):
        from .scene import vector
        p=vector(position)
        material=self.scene.material(self.scene.region_at(p))
        if material.kind!='dielectric' or self.scene.starts_terminated(p): return np.zeros(3,complex)
        return sum((s.field_at(p) for s in self.segments if s.material==material),start=np.zeros(3,complex))

    def fields_at(self, positions):
        return np.asarray([self.field_at(p) for p in positions])

    @property
    def power_balance(self):
        result = {"launched":self.launched_power,"absorbed":self.absorbed_power,
                  "interface_flux_residual":self.interface_flux_residual}
        for t in self.terminals: result[t.reason] = result.get(t.reason,0.0)+t.power
        return result


def interact(state,event,scene,frequency):
    """Fresnel splitting and differential (position, direction) transport.

    Curvature is geometric slope per transverse displacement. The interface
    map transforms both transverse axes, including oblique astigmatism.
    """
    d,u = state.direction,state.basis
    normal = event.normal.copy()
    if normal@d>0: normal = -normal
    c = -float(normal@d)
    if c < 1e-10: return [], "grazing"
    incoming = scene.material(state.region)
    obj = scene._objects[event.object_id]
    outgoing = obj.material if not obj.volume else scene.material(event.after)
    n1 = incoming.index(frequency).real
    n2 = outgoing.index(frequency).real if outgoing.kind != "pec" else 1.0
    rs,rp,ts,tp = fresnel(incoming,outgoing,c,frequency)
    s = np.cross(d,normal)
    s = unit(s) if np.linalg.norm(s)>1e-10 else u[:,0]
    pi = np.cross(s,d)
    es,ep = s@state.field, pi@state.field
    dn_dx = curvature(event.shape,event.position,normal)
    hit_projection = np.eye(3)-np.outer(d,normal)/(normal@d)
    dx = hit_projection@u
    dn = dn_dx@dx
    children = []
    for kind in ("reflect","transmit"):
        if kind == "transmit" and outgoing.kind != "dielectric": continue
        if kind == "reflect":
            direction = d+2*c*normal
            coefficient_s,coefficient_p = rs,rp
            derivative_x = -2*(np.outer(normal,d@dn)+(d@normal)*dn)
            derivative_d = (np.eye(3)-2*np.outer(normal,normal))@u
            region = state.region
        else:
            eta = n1/n2
            ct2 = 1-eta*eta*(1-c*c)
            if ct2 <= 1e-14:
                # Lossy evanescent transport near/above the critical angle is
                # outside the real-axis paraxial model; never invent a ray.
                if outgoing.index(frequency).imag or incoming.index(frequency).imag:
                    return [], "unsupported_lossy_critical"
                continue
            ct = math.sqrt(ct2)
            direction = eta*d+(eta*c-ct)*normal
            coefficient_s,coefficient_p = ts,tp
            dc = -d@dn
            derivative_x = (eta*c-ct)*dn+np.outer(normal,(eta-eta*eta*c/ct)*dc)
            derivative_d = eta*u-np.outer(normal,(eta-eta*eta*c/ct)*(normal@u))
            region = event.after
        direction = unit(direction)
        v = frame(direction)
        a = v.T@dx
        h = (v.T@derivative_x + v.T@derivative_d@state.curvature)@np.linalg.inv(a)
        if np.linalg.norm(h-h.T) > 1e-7*max(1,np.linalg.norm(h)):
            raise RuntimeError("Interface curvature lost reciprocity")
        h = (h+h.T)/2
        new_field = coefficient_s*es*s+coefficient_p*ep*np.cross(s,direction)
        ratio = float(np.vdot(new_field,new_field).real/max(np.vdot(state.field,state.field).real,1e-300))
        if kind == "transmit":
            adm1 = (incoming.index(frequency)/incoming.mu_r).real
            adm2 = (outgoing.index(frequency)/outgoing.mu_r).real
            ratio *= adm2*abs(direction@normal)/(adm1*c)
        children.append(State(event.position.copy(),direction,v,h,new_field,region,state.power*ratio,
                              state.distance,state.optical_path,
                              state.history+((event.object_id,kind,state.region,region),)))
    return children, ""


def solve(scene, source, config=None):
    """Trace one source or a coherent list of sources; sample via result.field_at.

    Finite beams use local paraxial interface maps. Aperture clipping and edge
    diffraction are not represented by a whole-beam split; warnings identify
    intersections whose footprints approach primitive edges.
    """
    from copy import copy
    cfg = copy(config) if config is not None else SolverConfig()
    cfg.validate()
    scene.validate()
    # Snapshot geometry so later caller edits cannot change saved field queries.
    snapshot = type(scene)(scene.background,tolerance=scene.tolerance)
    for o in scene.objects: snapshot._add(o.shape,o.material,o.volume,o.terminate)
    scene = snapshot
    for m in [scene.background]+[o.material for o in scene._objects if o.material and o.material.kind=="dielectric"]:
        n = m.index(cfg.frequency_hz)
        if n.imag/n.real > cfg.max_loss_ratio:
            raise ValueError("Dielectric loss exceeds configured weak-loss limit; strongly absorbing refraction is unsupported")
    sources = source if isinstance(source,(tuple,list)) else [source]
    if not sources: raise ValueError("At least one source is required")
    result = Simulation(cfg,scene)
    queue = []
    for src in sources:
        launched = src.launch(scene,cfg.frequency_hz)
        for state in launched if isinstance(launched,list) else [launched]:
            state.branch = len(queue)
            queue.append(state)
            result.launched_power += state.power
    if len(queue)>1:
        result.warnings.append("Coherent source superposition: power_balance sums constituent branch powers; use a field flux integral for physical total power")
    if len(queue)>cfg.max_branches: raise ValueError("Source count exceeds branch budget")
    allocated = len(queue)
    next_id = allocated
    def finish(state, reason):
        result.terminals.append(Terminal(state.branch,reason,state.position.copy(),state.power,state.history))
    while queue:
        state = queue.pop()
        if scene.starts_terminated(state.position):
            finish(state,"terminated"); continue
        if state.power == 0 or state.power < cfg.min_power_fraction*result.launched_power:
            finish(state,"power_cutoff"); continue
        remaining = cfg.max_distance-state.distance
        if remaining <= scene.tolerance:
            finish(state,"distance_limit"); continue
        event = scene.next_event(state.position,state.direction,remaining)
        length = remaining if event is None else event.distance
        material = scene.material(state.region)
        segment = Segment(state,length,material,cfg.frequency_hz)
        result.segments.append(segment)
        end = propagate(state,length,material,cfg.frequency_hz)
        result.absorbed_power += state.power-end.power
        if event is None:
            finish(end,"distance_limit"); continue
        if event.kind in ("ambiguous","terminate"):
            finish(end,"ambiguous" if event.kind=="ambiguous" else "terminated"); continue
        if len(state.history)>=cfg.max_interactions:
            finish(end,"interaction_limit"); continue
        # Matched-material boundaries may have changed region identity without
        # changing the field. Record the actual region at this interaction.
        end.region=event.before
        widths = segment.radii(length)
        if np.isfinite(widths).all():
            c = abs(state.direction@event.normal)
            footprint = 3*max(widths)/max(c,1e-10)
            if edge_distance(event.shape,event.position) < footprint:
                result.warnings.append(f"Branch {state.branch}: beam footprint approaches object {event.object_id} edge; clipping is unresolved")
        children,reason = interact(end,event,scene,cfg.frequency_hz)
        if reason:
            finish(end,reason); continue
        residual = end.power-sum(ch.power for ch in children)
        target = scene._objects[event.object_id]
        boundary_material = scene.material(event.after) if target.volume else target.material
        if boundary_material.kind in ('metal','pec'):
            result.absorbed_power += residual
        else:
            # In lossy incident media the interference flux of incident and
            # reflected waves prevents adding independent branch powers.
            result.interface_flux_residual += residual
            if abs(residual) > 1e-7*max(end.power,1e-300):
                result.warnings.append("Lossy-interface branch budgets omit interference flux; inspect interface_flux_residual")
        for child in children:
            child.parent,child.branch = state.branch,next_id
            next_id += 1
            if allocated >= cfg.max_branches:
                finish(child,"branch_limit")
            else:
                allocated += 1
                queue.append(child)
    result.warnings = list(dict.fromkeys(result.warnings))
    return result
