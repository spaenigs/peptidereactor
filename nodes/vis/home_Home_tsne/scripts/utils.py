from more_itertools import intersperse


def wrap_text(text, line_length=15):
    word = text.strip()
    words = word.split(" ")
    if len(words) > 1:
        wrapped_sentence = ""
        tmp_sentence = ""
        for w in words:
            if len(tmp_sentence + w) <= line_length:
                tmp_sentence += w + " "
            else:
                wrapped_sentence += tmp_sentence + "\n"
                tmp_sentence = w + " "
        return wrapped_sentence + tmp_sentence
    else:
        return "".join(intersperse("\n", word, n=line_length))