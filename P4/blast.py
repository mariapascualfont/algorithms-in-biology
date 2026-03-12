from Bio import SeqIO
from Bio.Align import substitution_matrices
import sys

subst_mat = substitution_matrices.load("BLOSUM62")

def build_index(db_seqs, k):
  db = {}
  for j in range(len(db_seqs)):
    for i in range(len(db_seqs[j])-k+1):
      if db_seqs[j][i:i+k] not in db:
        db[db_seqs[j][i:i+k]] = []
      db[db_seqs[j][i:i+k]].append((j,i))

  return dict(sorted(db.items()))

def find_seeds(query, db, k, t):

    query_kmers = [(i, query[i:i+k]) for i in range(len(query) - k + 1)]
    seeds = []
    for q_pos, q_word in query_kmers:
        for db_word, positions in db.items():
            joined_letters = zip(q_word, db_word)
            score = sum(subst_mat[a, b] for a, b in joined_letters)
            if score >= t:
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

def extend_seeds(seeds, query, db_seqs, k, X=15):

    hits = []  # will contain actual hits
    seen = (set())

    # to avoid duplicates, stores tuples (db_iseq, q_start, q_end, db_start, db_end)

    for s in seeds:
        q_pos = s["q_pos"]        # starting position in query
        db_iseq = s["db_iseq"]    # which db sequence
        db_pos = s["db_pos"]      # starting position in db sequence
        db_seq = db_seqs[db_iseq] # the actual db sequence string

        # --- Extend right ---
        score, best_score, best_right = 0, 0, 0
        right = 0
        while (
            q_pos + k + right < len(query)
            and db_pos + k + right < len(db_seq)
            and score >= best_score - X
        ):

            score += subst_mat[query[q_pos + k + right], db_seq[db_pos + k + right]]
            if score > best_score:
                best_score = score
                best_right = right + 1
            right += 1

        # --- Extend left ---
        score, best_score, best_left = 0, 0, 0
        left = 0
        while (
            q_pos - left - 1 >= 0 and db_pos - left - 1 >= 0 and score >= best_score - X
        ):
            score += subst_mat[query[q_pos - left - 1], db_seq[db_pos - left - 1]]
            if score > best_score:
                best_score = score
                best_left = left + 1
            left += 1

        q_start = q_pos - best_left
        q_end = q_pos + k + best_right - 1
        db_start = db_pos - best_left
        db_end = db_pos + k + best_right - 1

        key = (db_iseq, q_start, q_end, db_start, db_end)
        if key not in seen:
            seen.add(key)
            hits.append(
                {
                    "db_iseq": db_iseq,
                    "query_start": q_start,
                    "query_end": q_end,
                    "db_start": db_start,
                    "db_end": db_end,
                    "query_seq": query[q_start : q_end + 1],
                    "db_seq": db_seq[db_start : db_end + 1],
                    "seed": s["q_kmer"],
                }
            )

    return hits

def merge_overlapping(hits, query, db_seqs):

    merged_hits = []

    groups = {}
    for h in hits:
        diagonal = h["db_start"] - h["query_start"]  # we define the diagonal this way
        key = (h["db_iseq"], diagonal)
        if key not in groups:
            groups[key] = []
        groups[key].append(h)

    for (
        db_iseq,
        _diagonal,
    ), group in groups.items():  # iterate over each db sequence and diagonal
        group.sort(
            key=lambda x: x["query_start"]
        )  # sort hits by query_start to ensure they are in the correct order for merging
        current = dict(
            group[0]
        )  # start with the first hit in the group as the current hit to merge with subsequent hits

        for h in group[1:]:
            if h["query_start"] <= current["query_end"] + 1:
                current["query_end"] = max(current["query_end"], h["query_end"])
                current["db_end"] = max(current["db_end"], h["db_end"])  # similarly, the maximum end is selected from the hits
                current["query_seq"] = query[
                    current["query_start"] : current["query_end"] + 1
                ]
                current["db_seq"] = db_seqs[db_iseq][current["db_start"] : current["db_end"] + 1]
            else:
                merged_hits.append(current) # store the hit
                current = dict(h) # start a new current hit with the next hit in the group

        merged_hits.append(current)

    return merged_hits


if __name__ == "__main__":

    k = 3  # k-mer length
    t = 11  # BLOSUM62 score threshold

    query_file = r"C:\Users\maria\OneDrive\Escriptori\Universitat\Semester4\AB\P4\query.fasta"
    db_file =   r"C:\Users\maria\OneDrive\Escriptori\Universitat\Semester4\AB\P4\database.fasta"

    query = str(SeqIO.read(query_file, "fasta").seq)  # complete
    db_seqs = [str(record.seq) for record in SeqIO.parse(db_file, "fasta")]

    print(f"\nPreprocessing database (k={k})...")
    index = build_index(db_seqs, k)

    print(f"\nFinding seeds (threshold t={t})...")
    seeds = find_seeds(query, index, k, t)

    print(f"\nExtending seeds...")
    hits = extend_seeds(seeds, query, db_seqs, k)

    print(f"\nMerging overlapping hits...")
    merged = merge_overlapping(hits, query, db_seqs)

    print(f"\nMerged hits sorted by length (longest first): {len(merged)}")
    merged.sort(key=lambda h: h["query_end"] - h["query_start"], reverse=True)
    for h in merged[:10]:  # print top 10 longest hits
        diagonal = h["db_start"] - h["query_start"]
        print(
            f"  DB[{h['db_iseq']}] diag={diagonal}"
            f" Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_seq']}")
