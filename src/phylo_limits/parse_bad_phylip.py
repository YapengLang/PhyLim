# this file is from Gavin for parsing bad phylip files
import re
from collections import defaultdict
from cogent3.app.composable import define_app, LOADER
from cogent3.app import typing as c3types
from cogent3 import make_aligned_seqs, open_data_store

_wspace = re.compile(r"\s+")


@define_app(app_type=LOADER)
def load_bad_phylip(
    path: c3types.IdentifierType, moltype="dna"
) -> c3types.AlignedSeqsType:
    data = path.read()  # interim
    data = data.splitlines()
    num_seqs, seq_len = [int(e) for e in data[0].split()]
    labels = []
    seqs = defaultdict(list)
    seq_num = 0
    getting_labels = True
    for line in data[1:]:
        line = line.strip()
        if not line:
            continue

        if getting_labels:
            label, line = _wspace.split(line, maxsplit=1)
            labels.append(label)
            if len(labels) == num_seqs:
                getting_labels = False
        else:
            label = labels[seq_num]
            seq_num += 1
            if seq_num == num_seqs:
                seq_num = 0

        seq = _wspace.sub("", line)
        seqs[label].append(seq)

    seqs = {l: "".join(s) for l, s in seqs.items()}
    lengths = {len(s) for s in seqs.values()}
    assert lengths == {seq_len}
    return make_aligned_seqs(data=seqs, moltype=moltype)


if __name__ == "__main__":
    # works on Gavin's machine!
    dstore = open_data_store("~/Desktop/Inbox", suffix="phy", mode="r")
    loader = load_bad_phylip()
    # just doing one file
    print(repr(loader(dstore[0])))
