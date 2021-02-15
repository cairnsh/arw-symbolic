# HAC, 2020-10-15

from sympy import *
from progressbar import progressbar

def generate_progress_callback():
    pb = progressbar()
    def cb(what, col, *others):
        if what == "done":
            pb.display("done")
            pb.done()
            return
        out = "%s: column %d" % (what, col)
        if others:
            out += " " + str(others)
        pb.display(out)
    return cb

def inv(a, b, progress_callback=None, singular_callback=None):
    """Solve xa = b by doing column operations.
    This function destroys the original matrix a.
    
    It calls progress_callback as it goes:
      progress_callback("forward", i)
      progress_callback("backward", i)
      progress_callback("done", i)
    If progress_callback is None, this is displayed by a progressbar() so you can see how things are going.
    
    If a is singular, then:
      if singular_callback is None, the function dies with an Exception()
      otherwise,
        it calls singular_callback with the row where the pivoting failed and the partially lower-triangular matrix.
        This behaviour is used by nonzero_null_vector to extract a null vector.
    """
    x = b
    assert x.cols == a.cols
    if not progress_callback:
        progress_callback = generate_progress_callback()

    def SCALE(column, by):
        "Scale (column) by (by) and reduce fractions."
        a[:, column] = (by * a[:, column]).applyfunc(cancel)
        x[:, column] = (by * x[:, column]).applyfunc(cancel)

    def ADD(column_from, column_to, factor):
        if factor == 0:
            return
        "Add (by) * (column_from) to (column_to) and reduce fractions."
        a[:, column_to] = (a[:, column_to] + a[:, column_from] * factor).applyfunc(cancel)
        x[:, column_to] = (x[:, column_to] + x[:, column_from] * factor).applyfunc(cancel)

    def SWAP(column1, column2):
        "Swap (column1) and (column2)."
        if column1 == column2:
            return
        a[:, column1], a[:, column2] = a[:, column2], a[:, column1]
        x[:, column1], x[:, column2] = x[:, column2], x[:, column1]

    def DIE(i, error):
        print("failure in the inverter")
        print("debug information")
        print("i =", i)
        print("a[i, :] =", a[i, :])
        print("a = ", i)
        raise Exception(error)

    "Use column operations to make A lower-triangular with pivoting."
    "The pivot element is whichever one has the lowest degree."
    for i in range(a.cols):
        ip = None
        best_degree = None
        for ipivot in range(i, a.cols):
            if a[i, ipivot] != 0:
                d = total_degree(numer(a[i, ipivot]))
                if best_degree is None or d < best_degree:
                    ip = ipivot
                    best_degree = d
        if ip is None:
            # all of the entries a[i, j] for j >= i are zero, so we can't pivot.
            # this happens at some row if and only if a is singular
            if singular_callback:
                singular_callback(i, a)
                progress_callback and progress_callback("done", i)
                return
            else:
                DIE(i, "singular matrix")
        else:
            SWAP(i, ip) # move the pivot to the diagonal

        SCALE(i, 1/a[i, i])
        for j in range(i + 1, a.cols):
            ADD(i, j, -a[i, j])
            if progress_callback:
                progress_callback("forward", i, "adding", i, "to", j, "degree of last pivot:", best_degree)
        if progress_callback:
            progress_callback("forward", i, "degree of last pivot:", best_degree)

    "A is now upper-triangular; solve for x"
    for i in range(a.cols-1, -1, -1):
        for j in range(i):
            x[:, j] = (x[:,j] - x[:,i]*a[i,j]).applyfunc(cancel)
            if progress_callback:
                progress_callback("backward", i, "adding", i, "to", j)
            #for l in range(x.rows):
            #    x[l, j] = cancel(x[l, j] - x[l, i] * a[i, j])
        if progress_callback:
            progress_callback("backward", i)

    if progress_callback:
        progress_callback("done", -1)
    return x

def inverse(A, indices_from, indices_to):
    """Returns A^{-1}[indices_from, indices_to] as SparseMatrix.
    This function destroys the original matrix A."""
    assert A.rows == A.cols
    size = A.rows

    select = lambda what: {(j, what[j]): 1 for j in range(len(what))}
    b = SparseMatrix(len(indices_from), size, select(indices_from))
    x = inv(A, b)

    t = zeros(len(indices_from), len(indices_to))
    for i in range(len(indices_from)):
        for j in range(len(indices_to)):
            t[i, j] = x[i, indices_to[j]]
    return t

def nonzero_null_vector(A, progress_callback=None):
    """If A is singular, this function returns a nonzero right null vector.
    This function destroys the original matrix A."""
    null = [None]

    def singular_callback(i, a):
        # find a linear combination of the rows of a[0:i+1,0:i+1] that is zero
        coeffs = {i: 1}
        for j in range(i-1, -1, -1):
            coeffs[j] = cancel(-a[i, j]/a[j, j])
            # add -a[i, j]/a[j, j] times a[j, :] to a[i, :]
            for l in range(0, j+1):
                a[i, l] = cancel(a[i, l] + a[j, l] * coeffs[j])

        null[0] = SparseMatrix(1, size, {(0, j): coeffs[j] for j in range(i+1)})

    zero = SparseMatrix(0, size, {})
    inv(A, zero, progress_callback=progress_callback, singular_callback=singular_callback)
    return null[0]
