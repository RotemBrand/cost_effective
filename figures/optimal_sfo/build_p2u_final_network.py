from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.p2u_final_network import (
    BACKBONE_SUMMARY,
    FINAL_NETWORK_GPKG,
    FINAL_NETWORK_METADATA,
    build_and_write_final_network,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the final P2U optimized network after the ILP backbone.")
    parser.add_argument("--backbone-summary", type=Path, default=BACKBONE_SUMMARY)
    parser.add_argument("--output-gpkg", type=Path, default=FINAL_NETWORK_GPKG)
    parser.add_argument("--metadata-json", type=Path, default=FINAL_NETWORK_METADATA)
    parser.add_argument("--tree-mode", choices=["street_forest", "star"], default="street_forest")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = build_and_write_final_network(
        backbone_summary=args.backbone_summary,
        output_gpkg=args.output_gpkg,
        metadata_json=args.metadata_json,
        tree_mode=args.tree_mode,
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
