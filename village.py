


class Villages(object):
    """ A collection of Villages
    """

    def __init__(self, villages):
        """
        Parameters
        ----------
        villages : a list of figs.Village 
        """
        self.villages = villages
        self.current = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current >= len(self.villages):
            raise StopIteration
        else:
            self.current



class Village(object):
    """ Basic parameters for the village
    """

    def __init__(self, hostpopsize = 100,
            pervalence = 0.1):
        """
        """
        self.hostpopsize = 100
        self.prevalence = prevalence
    
