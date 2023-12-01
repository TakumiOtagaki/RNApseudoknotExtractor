# Base Pair Matrix (0 or 1) を受け取り、シュードノット（PK）を検知して渡す。
# 一つの塩基が複数の異なる箇所と塩基対を形成することは、無いと仮定している。


# this function can compress BasePairList.
# the output is compressed_BasePairList and the inv_hash table.
import argparse

def BasePairList_compression(BPL):
    # initialization:
    compressed_BPL = []
    inv_hash = [-1 for _ in range(2 * len(BPL))]
    hash = dict()
    # inv_hash[i] = < 小さい方から i 番目の塩基の実際のSeq での位置 (0-indexed)>
    # recursion:
    BaseSet = set()
    for (base1, base2) in BPL:
        BaseSet.update({base1, base2})
    sorted_Base = sorted(list(BaseSet))
    for i in range(len(sorted_Base)):
        inv_hash[i] = sorted_Base[i]
        hash[sorted_Base[i]] = i
    for (base1, base2) in BPL:
        compressed_BPL.append(
            (hash[base1], hash[base2]))
    num_of_bases_in_BPL = len(BaseSet)
    return compressed_BPL, inv_hash, num_of_bases_in_BPL

# df extraction.
# extract df["Saenger Classification"] == 19, 20, 28


def decompress_PKlayer_BPL(PK_layer, BPL, inv_hash):
    decompressed_PK_layer = []
    decompressed_BPL = []
    for (hashed_base1, hashed_base2) in PK_layer:
        decompressed_PK_layer.append(
            (inv_hash[hashed_base1], inv_hash[hashed_base2]))
    for (hashed_base1, hashed_base2) in BPL:
        decompressed_BPL.append(
            (inv_hash[hashed_base1], inv_hash[hashed_base2]))
    return decompressed_PK_layer, decompressed_BPL


def PKextractor(BPL, compression=True):
    PK_layers = []
    while BPL:
        # initialization:
        if compression:
            BPL, inv_hash, L = BasePairList_compression(BPL)
        else:
            L = max(BPL, key=lambda x: x[1])[1] + 1
        gamma = [[-1 for j in range(L)] for i in range(L)]
        for i in range(L):
            gamma[i][i] = 0
        for i in range(1, L):
            gamma[i][i-1] = 0

        # recursion:
        for d in range(1, L):  # k: diagonal index
            j = L
            for i in range(L-d-1, -1, -1):
                j -= 1
                max_candidate = [gamma[i+1][j], gamma[i][j-1],
                                 (gamma[i+1][j-1] + int((i, j) in BPL))]
                max_candidate.append(
                    max([(gamma[i][k] + gamma[k+1][j]) for k in range(i, j)]))
                gamma[i][j] = max(max_candidate)
        PK_layer, BPL = PK_traceback(gamma, BPL, L)
        if compression:
            PK_layer, BPL = decompress_PKlayer_BPL(PK_layer, BPL, inv_hash)
        PK_layers.append(PK_layer)
    return PK_layers


def PK_traceback(gamma, BPL, L):
    # initialization:
    route = ["." for _ in range(L)]
    trace_stack = [(0, L-1)]
    PK_layer = []  # base pair (i, j) を収納する。

    # recursion:
    while True:
        if not trace_stack:
            break
        i, j = trace_stack.pop(-1)
        if (i >= j):
            continue
        elif gamma[i+1][j] == gamma[i][j]:
            trace_stack.append((i+1, j))
        elif gamma[i][j-1] == gamma[i][j]:
            trace_stack.append((i, j-1))
        elif i+1 < L and gamma[i+1][j-1] + ((i, j) in BPL) == gamma[i][j]:
            if route[i] == "." and route[j] == ".":
                route[i], route[j] = "(", ")"
                basepair_ij = (i, j)
                PK_layer.append(basepair_ij)
                BPL.remove(basepair_ij)
                trace_stack.append((i+1, j-1))
        else:
            for k in range(i, j):
                if gamma[i][k] + gamma[k+1][j] == gamma[i][j]:
                    trace_stack += [(k+1, j), (i, k)]
                    break
    return PK_layer, BPL

def parser():
    parser = argparse.ArgumentParser(description="PKextractor")
    parser.add_argument("--input_file", type=str, help="input file path. 2 columns: i, j. Extension should be tsv file.")
    parser.add_argument("--output_file", type=str, help="output file path. 3 columns: i, j, layer_id. Extension should be tsv file.")
    args = parser.parse_args()
    return args.input_file, args.output_file


def main():
    input_file, output_file = parser()
    BPL = []
    with open(input_file, "r") as f:
        # input_file: tsv file. 2 columns. i,j means base pair (i,j).
        for line in f:
            i, j = line.split()
            BPL.append((int(i), int(j)))
    PK_layers = PKextractor(BPL)
    with open(output_file, "w") as f:
        # i, j, layer_id
        for PK_layer in PK_layers:
            for i, j in PK_layer:
                f.write(f"{i}\t{j}\t{PK_layers.index(PK_layer)}\n")


if __name__ == "__main__":
    main()

