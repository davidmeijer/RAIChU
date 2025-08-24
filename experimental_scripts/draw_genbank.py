import argparse 
import logging
import os


# install branch of raichu

# conda install "setuptools<81"
# pip install "setuptools<81"
import warnings
warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated",
    category=UserWarning
)

from raichu.representations import ClusterRepresentation
from raichu.run_raichu import build_cluster
from raichu.antismash import get_nrps_pks_modules


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True, help="out dir")
    return parser.parse_args()


def main():
    args = cli()
    logger.info(f"Input file: {args.i}")

    os.makedirs(args.o, exist_ok=True)

    module_blocks = get_nrps_pks_modules(args.i)

    cluster_representation = module_blocks.make_raichu_cluster()
    cluster_representation.write_cluster(args.o)
    cluster_representation = ClusterRepresentation.from_file(args.o)

    cluster = build_cluster(cluster_representation, strict=False)
    has_functional_modules = False
    for module in cluster.modules:
        if not module.is_broken:
            has_functional_modules = True
    if not has_functional_modules:
        exit(0)  # no functional modules, nothing to draw
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()

    print(cluster)
    for module in cluster.modules:
        print(module.type, module.subtype, module.is_starter_module, module.is_termination_module, module.tailoring_domains)
        # filter out relevant domains from (true, false) to false, true) and in between only (false, false)
        # give domains active/inactive label (see bgc 55)
        # just need to get coordinates gene in there so I can map sequence for parasect prediction
        # extract stereochem info like PKS_CIS and epimerization domains for the amino acids

    path = os.path.join(args.o, "cluster.svg")
    svg_str = cluster.draw_cluster(as_string=False, out_file=path)


if __name__ == "__main__":
    main()
