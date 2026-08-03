"""Small dependency-free validator for the contract schema features used here."""

from __future__ import annotations

import re
from typing import Any


class SchemaViolation(ValueError):
    pass


def validate(instance: Any, schema: dict[str, Any], root: dict[str, Any] | None = None, path: str = "document") -> None:
    root = root or schema
    if "$ref" in schema:
        target: Any = root
        for part in schema["$ref"].removeprefix("#/").split("/"):
            target = target[part]
        validate(instance, target, root, path); return
    if "anyOf" in schema:
        for option in schema["anyOf"]:
            try: validate(instance, option, root, path); return
            except SchemaViolation: pass
        raise SchemaViolation(f"{path}: value does not match any allowed form")
    if "const" in schema and instance != schema["const"]: raise SchemaViolation(f"{path}: must equal {schema['const']!r}")
    if "enum" in schema and instance not in schema["enum"]: raise SchemaViolation(f"{path}: unsupported value {instance!r}")
    types = schema.get("type")
    if types:
        types = [types] if isinstance(types, str) else types
        checks = {"object":lambda x:isinstance(x,dict),"array":lambda x:isinstance(x,list),"string":lambda x:isinstance(x,str),"integer":lambda x:isinstance(x,int) and not isinstance(x,bool),"number":lambda x:isinstance(x,(int,float)) and not isinstance(x,bool),"boolean":lambda x:isinstance(x,bool),"null":lambda x:x is None}
        if not any(checks[t](instance) for t in types): raise SchemaViolation(f"{path}: expected {' or '.join(types)}")
    if isinstance(instance, dict):
        for key in schema.get("required", []):
            if key not in instance: raise SchemaViolation(f"{path}: missing required property {key!r}")
        props=schema.get("properties", {})
        if schema.get("additionalProperties") is False:
            extra=set(instance)-set(props)
            if extra: raise SchemaViolation(f"{path}: unexpected property {sorted(extra)[0]!r}")
        for key,value in instance.items():
            child=props.get(key)
            if child is None and isinstance(schema.get("additionalProperties"),dict): child=schema["additionalProperties"]
            if child is not None: validate(value,child,root,f"{path}.{key}")
    if isinstance(instance,list):
        if len(instance)<schema.get("minItems",0): raise SchemaViolation(f"{path}: too few items")
        if "items" in schema:
            for index,value in enumerate(instance): validate(value,schema["items"],root,f"{path}.{index}")
    if isinstance(instance,str):
        if len(instance)<schema.get("minLength",0): raise SchemaViolation(f"{path}: string is too short")
        if "pattern" in schema and re.search(schema["pattern"],instance) is None: raise SchemaViolation(f"{path}: string does not match required pattern")
    if isinstance(instance,(int,float)) and not isinstance(instance,bool) and "minimum" in schema and instance<schema["minimum"]: raise SchemaViolation(f"{path}: below minimum")
