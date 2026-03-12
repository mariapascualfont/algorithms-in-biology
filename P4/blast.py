from Bio.Align import substitution_matrices

# from build_index import build_index

subst_mat = substitution_matrices.load("BLOSUM62")

def find_seeds(query, db, k, t):
    query_kmers = __ # complete
    seeds = []
    for __ # complete
        for __  # complete
            score = # complete
            if score # complete
                for db_iseq, db_pos in positions:
                    seeds.append(
                        {
                            "db_iseq": db_iseq,
                            "db_pos": db_pos,
                            "q_pos": q_pos,
                            "q_kmer": q_word,  # original query k-mer
                            "db_kmer": db_word,  # matched db k-mer (may differ)
                            "score": score,
                        }
                    )
    return seeds


if __name__ == "__main__":
    k = 3
    t = 11
    db = {
        "KTM": [(1, 0)],
        "MKT": [(1, 2)],
        "MRT": [(0, 0)],
        "RTA": [(0, 1)],
        "TAY": [(0, 2)],
        "TMK": [(1, 1)],
    }
    print(find_seeds("MKT", db, k, t))
