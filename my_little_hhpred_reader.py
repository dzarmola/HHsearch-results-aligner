import re

class HHpredHit:
    def __init__(self, query,hit_str):
        self.query = query
        self.target = ''
        self.prob = 0
        self.eval = 0
        self.score = 0
        self.aligned_cols = 0
        self.idents = 0
        self.q_cons = ''
        self.q_cons_s = 0 ## HHpred counts from 1
        self.q_cons_e = 0
        self.t_cons = ''
        self.t_cons_s = 0
        self.t_cons_e = 0

        self.read_in(hit_str)

        self.match()

    def __str__(self):
        return "{t}: E-value: {e} Q:{qs}-{qe} T:{ts}-{te}".format(t=self.target,
                        e=self.eval, qs=self.q_cons_s,qe=self.q_cons_e,ts=self.t_cons_s,
                        te=self.t_cons_e)

    def __repr__(self):
        return self.__str__()

    def read_in(self,txt):
        self.target = re.search(">([A-Za-z0-9_]+)",txt).group(1)
        self.prob = float(re.search("Probab=([0-9\.]+)",txt).group(1))
        self.eval = float(re.search("E-value=([0-9\.e+-]+)",txt).group(1))
        self.score = float(re.search("Score=([0-9\.]+)",txt).group(1))
        self.aligned_cols = int(re.search("Aligned_cols=([0-9\.]+)",txt).group(1))
        self.idents = float(re.search("Identities=([0-9\.]+)%",txt).group(1))
        q_cons = zip(*re.findall("Q Consensus +([0-9]+) (.*?) +([0-9]+)",txt))
        q_ref = zip(*re.findall("Q "+self.query[:14]+" +([0-9]+) (.*?) +([0-9]+)",txt))
        self.q_ref = "".join(q_ref[1])
        self.q_cons = "".join(q_cons[1])
        self.q_cons_s = int(q_cons[0][0])
        self.q_cons_e = int(q_cons[2][-1])
        t_cons = zip(*re.findall("T Consensus +([0-9]+) (.*?) +([0-9]+)",txt))
        t_ref = zip(*re.findall("T "+self.target[:14]+" +([0-9]+) (.*?) +([0-9]+)",txt))
        self.t_cons = "".join(t_cons[1])
        self.t_ref = "".join(t_ref[1])
        self.t_cons_s = int(t_cons[0][0])
        self.t_cons_e = int(t_cons[2][-1])

        self.q_better_ref = "".join([self.q_ref[x] for x in xrange(len(self.q_ref)) if self.q_cons[x] not in "-."])
        self.t_better_ref = "".join([self.t_ref[x] for x in xrange(len(self.t_ref)) if self.t_cons[x] not in "-."])

    def match(self):
        self.match = []
        ci = self.q_cons_s
        cj = self.t_cons_s
        for i,j in zip(self.q_cons, self.t_cons):
            p = [None,None]
            if i not in ["-","."]:
                p[0] = ci
                ci+=1
            if j not in ["-","."]:
                p[1] = cj
                cj+=1
            self.match.append(p)

class HHpredOutput:
    def __init__(self,filename):
        self.file = filename
        self.query = ''
        self.hits = []

        self.read_in_data()

    def read_in_data(self):
        with open(self.file) as input:
            self.query = re.search("Query\s+([A-Za-z0-9_]+)",input.readline()).group(1)
            self.len = int(re.search("Match_columns\s+([0-9]+)",input.readline()).group(1))
            _hits = re.split("No [0-9]+",input.read())
            for _h in _hits[1:]:
                self.hits.append(HHpredHit(self.query,_h))
                if self.hits[-1].target == self.query:
                    self.hits.pop()


if __name__ == "__main__":
    import sys
    h = HHpredOutput(sys.argv[1])
    print h.query
    print h.hits
    print h.hits[0].q_better_ref
    print h.hits[0].t_better_ref
    print h.hits[-1].match
