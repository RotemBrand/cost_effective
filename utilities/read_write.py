import networkx as nx
import pandas as pd
import numpy as np
import os

import json
import warnings
from pathlib import Path

def write_nxjson(df: pd.DataFrame, path: str | Path):
    # transform df
    df_to_save = df.map(_nx_to_json)
    df_to_save.reset_index(drop=True, inplace=True)

    # save
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df_to_save.to_json(path, orient='records', lines=True)
    print(f"SAVE DATA TO {path}")

def read_nxjson(path: str | Path) -> pd.DataFrame:
    df = pd.read_json(path, lines=True)
    df = df.map(_json_to_nx)
    return df


def _nx_to_json(graph):
    if graph.__class__.__name__ in NAME_TO_CLASS:
        return {
        'type': graph.__class__.__name__,
        'graph': {
            'nodes': [[n, data] for n, data in graph.nodes.items()],
            'edges': [[e, data] for e, data in graph.edges.items()],
            'graph': [[k, v] for k, v in graph.graph.items()]
        }
    }
    return graph

def _json_to_nx(obj):
    # validate if nxjson
    is_nx_json = isinstance(obj, dict) and (obj.get('type') in NAME_TO_CLASS)
    if not is_nx_json:
        return obj

    # construct graph
    g_type = NAME_TO_CLASS[obj['type']]
    graph =  g_type()
    for node, data in obj['graph']['nodes']:
        graph.add_node(node, **data)
    for edge, data in obj['graph']['edges']:
        graph.add_edge(*edge, **data)
    if 'graph' in obj['graph']:
        graph.graph = dict(obj['graph']['graph'])
    return graph


NAME_TO_CLASS = {
    'Graph' : nx.Graph,
    "DiGraph" : nx.DiGraph,
    'MultiGraph' : nx.MultiGraph,
    "MultiDiGraph" : nx.MultiDiGraph
}




################################ old ##############
def pd_to_json(df: pd.DataFrame, file_name: str):
    warnings.warn("dont use pd_to_json", category=DeprecationWarning, stacklevel=2)
    dict_to_json(df.reset_index().to_dict(), file_name)

def pd_read_json(file_name: str) -> pd.DataFrame:
    warnings.warn("dont use pd_read_json", category=DeprecationWarning, stacklevel=2)
    js = dict_read_json(file_name)
    df = pd.DataFrame(js)
    return df

def dict_to_json(d: dict, file_name: str):
    """
    Recursively convert dict values to JSON-serializable,
    detecting NetworkX graphs, and save to a file.
    """
    def value_to_json(v):
        if isinstance(v, dict):
            return {k: value_to_json(x) for k, x in v.items()}
        if isinstance(v, list):
            return [value_to_json(x) for x in v]
        if v.__class__.__name__ in NAME_TO_CLASS:
            return nx_to_json_str(v)
        return v

    js = value_to_json(d)
    s = json.dumps(js, ensure_ascii=False, indent=2)
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    with open(file_name, "w", encoding="utf-8") as f:
        f.write(s)

def dict_read_json(file_name: str):
    warnings.warn("dont use dict_read_json", category=DeprecationWarning, stacklevel=2)
    def value_to_nx(x):
        if _is_nx_dict(x):
            return value_to_nx(x)
        if isinstance(x, dict):
            return {convert_key(k): value_to_nx(v) for k, v in x.items()}
        if isinstance(x, list):
            return [value_to_nx(v) for v in x]
        if isinstance(x, str) and x.startswith('{'):
            try:
                return json_to_nx_old(json.loads(x))
            except:
                return x
        return x
    # read
    with open(file_name, "r", encoding="utf-8") as f:
        data = json.load(f)
    return value_to_nx(data)


def _is_nx_dict(x) -> bool:
    if isinstance(x, dict) and 'graph' in x and 'nodes' in x and 'edges' in x:
        return True
    return False

def convert_key(k: str):
    """
    Try to convert a JSON object key (always str) into
    int, float, bool, or None if possible. Otherwise return as str.
    """
    if not isinstance(k, str):
        return k

    # bool
    if k.lower() == "true":
        return True
    if k.lower() == "false":
        return False

    # None
    if k.lower() == "null" or k.lower() == "none":
        return None

    # int
    try:
        return int(k)
    except ValueError:
        pass

    # float
    try:
        return float(k)
    except ValueError:
        pass

    # fallback: keep as str
    return k

def safe_key(x):
    try:
        # --- numpy types ---
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, (np.integer,)):
            return int(x)
        if isinstance(x, (np.floating,)):
            return float(x)


        # --- already JSON-serializable? ---
        json.dumps(x)
        return x
    except (TypeError, OverflowError):
        # fallback: drop or convert to string
        return str(x)

def safe_attrs(data: dict) -> dict:
    """
    Convert a dict of attributes into something JSON-serializable.
    Handles common non-JSON types (numpy, sets, tuples).
    Skips anything else that cannot be converted.
    """
    out = {}
    for k, v in data.items():
        try:
            # --- numpy types ---
            if isinstance(v, np.ndarray):
                out[k] = v.tolist()
                continue
            if isinstance(v, (np.integer,)):
                out[k] = int(v)
                continue
            if isinstance(v, (np.floating,)):
                out[k] = float(v)
                continue

            # --- sets / tuples ---
            if isinstance(v, (set, tuple)):
                out[k] = list(v)
                continue

            # --- already JSON-serializable? ---
            json.dumps(v)
            out[k] = v
        except (TypeError, OverflowError):
            # fallback: drop or convert to string
            out[k] = str(v)
            continue
    return out

def nx_to_json_str(graph) -> str:
    return json.dumps({
        'type': graph.__class__.__name__,
        'graph': {
            'nodes': [[safe_key(n), safe_attrs(data)] for n, data in graph.nodes.items()],
            'edges': [[str(e), safe_attrs(data)] for e, data in graph.edges.items()],
            'graph': safe_attrs(graph.graph)
        }
    })


def json_to_nx_old(js) -> nx.Graph:
    g_type = NAME_TO_CLASS[js['type']]
    graph =  g_type()
    for node, data in js['graph']['nodes']:
        graph.add_node(node, **data)
    for edge, data in js['graph']['edges']:
        if isinstance(edge, str):
            edge = eval(edge)
        graph.add_edge(*edge, **data)
    if 'graph' in js['graph']:
        graph.graph = js['graph']['graph']
    return graph
