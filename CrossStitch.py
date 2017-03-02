from string import ascii_lowercase
import sys


def eprint(s="", end="\n"):
    sys.stderr.write(s + end)



class CrossStitch:
    def embroider(self, pattern):

        eprint(str(pattern))

        ret = []
        S = len(pattern)
        for col in ascii_lowercase:
            eprint(col)
            first = True
            for r in range(S):
                for c in range(S):
                    if pattern[r][c] == col:
                        if first:
                            first = False
                            ret.append(col)
                        ret.append(str(r+1) + " " + str(c))
                        ret.append(str(r) + " " + str(c+1))
                        ret.append(str(r+1) + " " + str(c+1))
                        ret.append(str(r) + " " + str(c))

        for r in ret:
            if r == "a" or r == "b":
                eprint("\n" + r)
            else:
                eprint(r, " ")
        eprint()
        eprint(str(ret))

        return ret

# -------8<------- end of solution submitted to the website -------8<-------

import sys
S = int(input())
pattern = []
for i in range(S):
    pattern.append(input().strip())

cs = CrossStitch()
ret = cs.embroider(pattern)
print(len(ret))
for st in ret:
    print(st)
sys.stdout.flush()




