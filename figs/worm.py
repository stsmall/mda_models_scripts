

class Worms(object):
    def __init__(self, meta, haplotype1=None, haplotype2=None):
        self.meta = meta 
        if haplotype1:
            self.h1= {}
        else:
            self.h1 = haplotype1
        if haplotype2:
            self.h2 = {}
        else:
            self.h2 = haplotype2

