import argparse 
import os

from raichu.retromol import (
    draw_cluster,
    gbk_to_cluster_files, 
    get_biosynthetic_modules,
    read_cluster_files, 
)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True, help="out dir")
    return parser.parse_args()


def main():
    args = cli()
    os.makedirs(args.o, exist_ok=True)

    cluster_file_core = os.path.join(args.o, "cluster.txt")
    gbk_to_cluster_files(args.i, cluster_file_core)
    cluster = read_cluster_files(cluster_file_core)
    if not cluster:
        exit(0)

    biosynthetic_modules = get_biosynthetic_modules(cluster)
    if biosynthetic_modules:
        for m in biosynthetic_modules:
            print(m.id, m.subtype, m.is_starter_module, m.is_termination_module)
            for d in m.domains:
                print(type(d))

    path = os.path.join(args.o, "cluster.svg")
    svg_str = draw_cluster(cluster)
    with open(path, "w") as f:
        f.write(svg_str)


if __name__ == "__main__":
    main()



# need in biosynthetic_modules:
# - stereochem pks modules, epimerization for amino acids or not
# - substrate (methylmalonyl or malonyl or ... and amino acid sequence for A-domains)
# - read out of biosynthetic type so I can filter out iterative PKS
# - put general information about organism at top of cluster file