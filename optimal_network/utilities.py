import time

TIME = 0
def toc(s: str='', debug: bool=True):
    if debug:
        global TIME
        if s != '':
            print(f'{s}: {time.time() - TIME}')
        TIME = time.time()


class ILPNoSolutionError(Exception):
    """Raised if not fassible solution found in an ILP problem"""
    pass

