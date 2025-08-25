import argparse 
import os

# import retromol module first because it suppresses deprecation warning setuptools
from raichu.retromol import (
    draw_cluster,
    gbk_to_cluster_files, 
    get_biosynthetic_modules,
    read_cluster_files, 
)
from raichu.module import ModuleType
from raichu.domain.domain_types import TailoringDomainType as TD, SynthesisDomainType as SD, RecognitionDomainType as RD


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True, help="out dir")
    return parser.parse_args()


def main():
    args = cli()
    os.makedirs(args.o, exist_ok=True)

    cluster_file_core = os.path.join(args.o, "cluster.txt")
    cluster_file_tailoring = os.path.join(args.o, "tailoring.txt")

    # TODO: get general organism information into cluster file, make sure it can be read as well

    gbk_to_cluster_files(args.i, cluster_file_core, cluster_file_tailoring)
    cluster = read_cluster_files(cluster_file_core)
    if not cluster:
        exit(0)
    
    # Experiment with readout
    biosynthetic_modules = get_biosynthetic_modules(cluster)
    # biosynthetic_modules = cluster.modules

    primary_sequence = []
    if biosynthetic_modules:
        for m in biosynthetic_modules:
            # print(m.type)
            # if m.type == ModuleType.PKS and any(d.type for d in m.domains if d.type == SD.KS and d.active and d.used):
            if m.type == ModuleType.PKS and any(d.type for d in m.domains if d.active and d.used):
                tailoring_domains = [d.type for d in m.domains if d.type in [TD.KR, TD.DH, TD.ER] and d.active and d.used]
                expected_domains_A = []
                expected_domains_B = [TD.KR]
                expected_domains_C = [TD.KR, TD.DH]
                expected_domains_D = [TD.KR, TD.DH, TD.ER]
                # if there is no KS it is a starter so "other"
                if not any(d.type == SD.KS for d in m.domains if d.active and d.used):
                    primary_sequence.append("other")
                    continue
                # if type sets are identical assign type A/B/C/D
                if set(tailoring_domains) == set(expected_domains_A):
                    motif_name = "A"
                elif set(tailoring_domains) == set(expected_domains_B):
                    motif_name = "B"
                elif set(tailoring_domains) == set(expected_domains_C):
                    motif_name = "C"
                elif set(tailoring_domains) == set(expected_domains_D):
                    motif_name = "D"
                else:
                    print(tailoring_domains)
                    motif_name = "A"

                at_domains = [d for d in m.domains if d.type == RD.AT]
                if len(at_domains) != 1:
                    primary_sequence.append(motif_name)
                    continue
                at_domain = at_domains[0]
                # WILDCARD/MALONYL_COA/METHYLMALONYL_COA/METHOXYMALONYL_COA/ETHYLMALONYL_COA
                # print(at_domain.type.name, at_domain.substrate.name, [d.name for d in tailoring_domains]) 

                substrate_name = at_domain.substrate.name
                if substrate_name == "WILDCARD":
                    pass
                elif substrate_name == "MALONYL_COA":
                    unit_name = f"{motif_name}1"
                elif substrate_name == "METHYLMALONYL_COA":
                    unit_name = f"{motif_name}2"
                elif substrate_name == "METHOXYMALONYL_COA":
                    unit_name = f"{motif_name}6"
                elif substrate_name == "ETHYLMALONYL_COA":
                    unit_name = f"{motif_name}4"
                else:
                    unit_name = motif_name
                
                for d in m.domains:
                    print(d.type, d.subtype)
                print(unit_name)
                primary_sequence.append(unit_name)

                # TODO: get info on stereochemistry
                # for KR: A1, A2, B1, B2, C1, C2
                #       A1/A2/B1/B2 different stereoisomers
                #       C1 inactive
                #       C2 no ketoreductase activity but does catalyze epimerization
                #       SEE RAICHU FIGURE S2 FOR IMPLEMENTATION
                # for ER: R, S, UNKNOWN

            elif m.type == ModuleType.NRPS:
                domain_types = [d.type for d in m.domains if d.active and d.used]
                if any(d == RD.A for d in domain_types):
                    has_epimerization = any(d.type == TD.E for d in m.domains if d.active and d.used)

                    a_domains = [d for d in m.domains if d.type == RD.A and d.active and d.used]
                    assert len(a_domains) == 1, "Expected exactly one A-domain"
                    a_domain = a_domains[0]
                    # print(type(a_domain))
                    # print(a_domain.translation)
                    primary_sequence.append("amino_acid")

                    # TODO: epimerization of α-carbon of l-amino acid substrates

            else:
                pass
    print(primary_sequence)

    path = os.path.join(args.o, "cluster.svg")
    svg_str = draw_cluster(cluster)
    with open(path, "w") as f:
        f.write(svg_str)


if __name__ == "__main__":
    main()

