"""Versioned, portable JSON results. No pickle and no legacy schema migration."""
import json
from pathlib import Path
from dataclasses import fields, is_dataclass
import numpy as np
from .materials import Material
from .scene import Scene, Object, Box, Rectangle, Disk, Sphere, Cylinder, Triangle
from .beams import State, Segment
from .transport import SolverConfig, Terminal, Simulation

_types = {t.__name__:t for t in (Material,Object,Box,State,Segment,SolverConfig,Terminal)}
_shape_fields = {
    Rectangle:("center","u","v","half_width","half_height"),
    Disk:("center","normal","radius"), Sphere:("center","radius"),
    Cylinder:("center","axis","radius","half_length","capped"), Triangle:("a","b","c")}


def _encode(value):
    if isinstance(value,np.ndarray): return {"type":"array","value":_encode(value.tolist())}
    if isinstance(value,complex): return {"type":"complex","real":value.real,"imag":value.imag}
    if isinstance(value,np.generic): return _encode(value.item())
    if type(value) in _shape_fields:
        return {"type":type(value).__name__,"args":[_encode(getattr(value,k)) for k in _shape_fields[type(value)]]}
    if is_dataclass(value):
        return {"type":type(value).__name__,"fields":{f.name:_encode(getattr(value,f.name)) for f in fields(value)}}
    if isinstance(value,(list,tuple)): return [_encode(v) for v in value]
    if isinstance(value,dict): return {k:_encode(v) for k,v in value.items()}
    return value


def _decode(value):
    if isinstance(value,list): return [_decode(v) for v in value]
    if not isinstance(value,dict): return value
    kind = value.get("type")
    if kind == "array": return np.asarray(_decode(value["value"]))
    if kind == "complex": return complex(value["real"],value["imag"])
    shapes = {t.__name__:t for t in _shape_fields}
    if kind in shapes: return shapes[kind](*[_decode(v) for v in value["args"]])
    if kind in _types:
        args = {k:_decode(v) for k,v in value["fields"].items()}
        if "history" in args: args["history"] = tuple(tuple(x) for x in args["history"])
        return _types[kind](**args)
    if kind is not None: raise ValueError(f"Unsupported result record {kind}")
    return {k:_decode(v) for k,v in value.items()}


def save_result(result, path):
    data = {"schema":"raytracing.beams","version":1,"phase":"exp(-i omega t)",
            "background":result.scene.background,"tolerance":result.scene.tolerance,
            "objects":result.scene.objects,"config":result.config,"segments":result.segments,
            "terminals":result.terminals,"warnings":result.warnings,
            "launched_power":result.launched_power,"absorbed_power":result.absorbed_power,
            "interface_flux_residual":result.interface_flux_residual}
    Path(path).write_text(json.dumps(_encode(data),indent=2,allow_nan=False),encoding="utf-8")


def load_result(path):
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    if raw.get("schema") != "raytracing.beams" or raw.get("version") != 1 or raw.get("phase") != "exp(-i omega t)":
        raise ValueError("Unsupported result schema or phasor convention")
    data = _decode(raw)
    scene = Scene(data["background"],tolerance=data["tolerance"])
    for o in data["objects"]: scene._add(o.shape,o.material,o.volume,o.terminate)
    data["config"].validate()
    return Simulation(data["config"],scene,data["segments"],data["terminals"],data["warnings"],
                      data["launched_power"],data["absorbed_power"],data["interface_flux_residual"])
