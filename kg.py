from collections import defaultdict
from itertools import chain
from operator import itemgetter
import logging
from pprint import pprint
import re

from plp import PLP


p = PLP()


def is_any_cz_mowy(word, cz_mowy_s):
    ks = p.rec(word)
    if not ks:
        return True
    return any(is_cz_mowy(ks, cz_mowy)
               for cz_mowy in cz_mowy_s)


def is_cz_mowy(ks, cz_mowy):
    return any(p.label(k)[0] == cz_mowy
               for k in ks)


def read_articles():
    inside_block = False
    block_start = re.compile(r'<doc id="\d+" url=".*" title=".*">')
    block_end = re.compile(r'</doc>')
    doc = []
    with open("../test_ss.txt", "r", encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            if i % 109558 == 0:
                logging.debug("Progress: {0:d}%".format(i // 109558))
            if not inside_block:
                if block_start.match(line):
                    inside_block = True
            elif block_end.match(line):
                yield ''.join(doc)
                doc = []
                inside_block = False
            else:
                doc.append(line)
        if doc:
            yield ''.join(doc)


def preparse_lines(inp):
    transformations = [(re.compile(r"'s"),  # joe's -> joe
                        ''),
                       (re.compile(r'/'),  # TCP/IP -> TCPIP
                        ''),
                       (re.compile(r'\b[0-9]+\b'),  # usuwanie liczb
                        ''),
                       (re.compile(r'[^a-zA-Z0-9żółćęśąźńŻÓŁĆĘŚĄŹŃ]+'),
                        ' ')]
    for doc in inp:
        for reg, repl in transformations:
            doc = reg.sub(repl, doc)
        doc = doc.lower()
        yield doc


def groups_filter(lst, pred):
    lst = list(lst)
    len_lst, i = len(lst), 0
    while i < len_lst:
        while i < len_lst and not pred(lst[i]):
            i += 1
        if i < len_lst:
            k = i
            while k < len_lst and pred(lst[k]):
                k += 1
            yield lst[i:k]
            i = k
        else:
            yield lst[i:]
            break


def words_ok(word):
    words_ok_roman = re.compile(r'^[ivx]+$')
    try:
        int(word)
    except ValueError:
        pass
    else:
        return False

    if words_ok_roman.match(word):
        return False

    return is_any_cz_mowy(word, [
        # PLP.CZESCI_MOWY.CZASOWNIK,
        # PLP.CZESCI_MOWY.LICZEBNIK,
        # PLP.CZESCI_MOWY.NIEODMIENNY,
        PLP.CZESCI_MOWY.PRZYMIOTNIK,
        # PLP.CZESCI_MOWY.PRZYSLOWEK,
        PLP.CZESCI_MOWY.RZECZOWNIK,
        # PLP.CZESCI_MOWY.SEGMENTOWY,
        # PLP.CZESCI_MOWY.SKROT,
        # PLP.CZESCI_MOWY.ZAIMEK
    ])


def group_forms(inp):
    for doc in inp:
        for x in groups_filter(doc.split(), words_ok):
            if len(x) > 1:
                yield x


def to_base(word):
    base = p.rec(word)
    if base:
        return p.bform(base[0])
    else:
        return word


def mk_dict(inp):
    res = defaultdict(lambda: defaultdict(int))
    for doc in inp:
        try:
            previous, tail = to_base(doc[0]), doc[1:]
            for current in tail:
                curr = to_base(current)
                res[previous][curr] += 1
                previous = curr
        except IndexError:
            pass
    return dict(res)


def calc_score(dct):
    for word, nexts in list(dct.items()):
        assert isinstance(word, str)
        assert isinstance(nexts, dict)
        if max(nexts.values()) < 3:
            del dct[word]
        else:
            divider = sum(nexts.values())
            dct[word] = {word: (score / divider)
                         for word, score
                         in nexts.items()
                         if score > 1}
    return dct


def mk_sorted(dct):
    dct_a = {}
    for word, nexts in dct.items():
        dct_a[word] = sorted(((score, [word])
                              for word, score in nexts.items()),
                             key=itemgetter(0),
                             reverse=True)
    return dct, dct_a


def baseline_filter_scored(dct_a):
    for previous, arr in dct_a.items():
        for score, current in arr[:len(arr) // 2]:
            yield (score, previous, current)


def get_sorted_lst(lst):
    return sorted(lst,
                  key=itemgetter(0),
                  reverse=True)


def generate_longer_seqs(dct, dct_a, markov_maxlen=5):
    for i in range(1, markov_maxlen - 1):
        for X_head, X_tails in dct_a.items():
            new_items = []
            for X_score, X_tail in X_tails:
                if (X_score > 0.3
                        and len(X_tail) == i
                        and X_tail[-1] in dct):
                    for Y, xY_score in dct[X_tail[-1]].items():
                        # noinspection PyPep8Naming
                        XY_score = (xY_score + X_score) * 0.5
                        if XY_score > 0.4:
                            print(">>", XY_score, X_tail + [Y])
                            new_items.append((XY_score, X_tail + [Y]))
            X_tails.extend(new_items)
    return dct_a


def trunc_dcta(dct_a):
    dct = {X_head: {tuple(Y): [score, defaultdict(int)]
                    for score, Y in X_successors
                    if score > 0.1}
           for X_head, X_successors in dct_a.items()}
    return dct


def inin_inc(elem_a, elem_b, ddct, sentence):
    try:
        ddct[elem_a][elem_b][1][sentence] += 1
    except KeyError:
        pass


def get_type(word):
    # noinspection PyBroadException
    try:
        return p.label(p.rec(word)[0])[0]
    except:
        return '?'


def get_statements(dct_p):
    rule = re.compile(r'^(AC*)+$')
    for X_head, X_succs in dct_p.items():
        for X_tail, X_tail_res in X_succs.items():
            # noinspection PyPep8Naming
            [score, X_tail_dct] = X_tail_res
            if score > 0.5 and X_tail_dct:
                for word, times in X_tail_dct.items():
                    if rule.match(''.join(get_type(x) for x in word.split())):
                        yield word, score * times


def main():
    ____log = logging

    stream = read_articles()
    stream = preparse_lines(stream)
    stream = group_forms(stream)
    ____log.info("Done making streams")

    dct = mk_dict(stream)
    ____log.info("Done making initial dict")

    dct = calc_score(dct)
    ____log.info("Done calculating scores")

    dct, dct_a = mk_sorted(dct)
    ____log.info("Done sorting results")

    dct_a = generate_longer_seqs(dct, dct_a)
    ____log.info("Done Markoving")

    dct_p = trunc_dcta(dct_a)  # :: Map str
                               #        (Map (Tuple str*)
                               #             [int,
                               #              (Map str int)])
    pprint(dct_p)
    stream = read_articles()
    stream = preparse_lines(stream)
    stream = map(lambda x: x.split(), stream)
    stream = chain.from_iterable(stream)
    a, b, c, d, e = '', '', '', '', ''
    ab, bb, cb, db, eb = '', '', '', '', ''
    for word in stream:
        a, b, c, d, e = b, c, d, e, word
        abcde = [a, b, c, d, e]
        ab, bb, cb, db, eb = bb, cb, db, eb, to_base(word)
        inin_inc(ab, (bb, cb, db, eb), dct_p, ' '.join(abcde))
        inin_inc(bb, (cb, db, eb), dct_p, ' '.join(abcde[1:]))
        inin_inc(cb, (db, eb), dct_p, ' '.join(abcde[2:]))
        inin_inc(db, (eb, ), dct_p, ' '.join(abcde[3:]))

    stream = filter(lambda x: x[1] > 1, get_statements(dct_p))
    print(sorted(stream, key=itemgetter(1), reverse=True))

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()
