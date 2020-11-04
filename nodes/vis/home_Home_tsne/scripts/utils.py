from functools import reduce
from more_itertools import intersperse


def wrap_text(text, line_length=15):
    def concat(old, new: str):
        if len(new) > line_length:
            old += [new]
        elif len(old) == 0 or \
                len(old[-1]) + len(new) >= line_length:
            old += [f"{new} "]
        else:
            old[-1] += new
        return old
    words = text.strip().split(" ")
    if len(words) == 1:
        return "".join(intersperse("\n", words[0], n=line_length))
    else:
        parts = reduce(lambda res, w: concat(res, w), words, [])
        return "\n".join(parts[:8])
