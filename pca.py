#!/usr/bin/env python3

import random, math
import matplotlib.pyplot as plot
import operator
import itertools

"""
    Principal Component Analysis
    Here:
    project 2d dataset onto a one dimension vector
"""

"""
    random point
"""
def random_point(mini, maxim):
    x = random.randrange(mini, maxim)
    y = random.randrange(mini, maxim)
    return (x,y)

"""
    fill the list with non-overlapping random points
    pool (hashtable) is used to check for overlapping
"""
def fill(l, n, mini, maxim):
    length = n
    point_pool = {}
    attribute = {}
    mean_x, mean_y = 0,0
    while n>0:
        (x,y) = random_point(mini, maxim)
        already = point_pool.get( (x,y), False)
        if not already:
            mean_x += x
            mean_y += y
            l.append( (x,y) )
            n -= 1
            point_pool[ (x,y) ] = True
        else:
            point_pool[ (x,y) ] = False
            continue
    del point_pool
    l.append( (mean_x/length, mean_y/length) )
    return (mean_x/length, mean_y/length)

"""
    mean normalization : 
    shift mean to origin
"""
def normalize(l, mean):
    sum_var = 0
    for index, point in enumerate(l):
        l[index] = (l[index][0] - mean[0], l[index][1]-mean[1])
        #l[index] = tuple(map(lambda x, y: x - y, l[index], mean))

"""
    find mean
"""
def mean_calc(l):
    n = len(l)
    x,y = zip(*l)
    return ( sum(x)/n, sum(y)/n)

""" 
    find covariance of data set x,y
"""
def covariance(xl, yl, mean=(0,0) ):
    n = len(xl)
    summation = 0
    for i in range(n):
        summation += (xl[i]-mean[0]) * (yl[i]-mean[1])
    #print (sum( [ (x-mean[0]) * (y-mean[1]) for x,y in zip(xl, yl)] ) / (n-1))
    return summation / (n-1)

"""
    calculate the covariance matrix
    mean is pre calculated during normalization or filling
"""
def covariance_matrix(l, mean):
    flattened = list(zip(*l))
    n = len(flattened)
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            value = covariance(flattened[i], flattened[j], mean)
            row.append(value)
        matrix.append(row)
    """
    matrix = [
                [ covariance(xl,xl,mean), covariance(xl, yl, mean)],
                [ covariance(yl,xl,mean), covariance(yl, yl, mean)],
            ]
    """
    return matrix

"""
    matrix * vector
"""
def multiply_eigen(mat, vec, size):
    res = [1] * size
    for i in range(size):
        sum = 0
        for j in range(size):
            sum += mat[i][j] * vec[j]
        res[i] = sum
    great = abs(max(res))
    res = [ (x/great) for x in res ]
    return res, great

def multiply(mat, vec, size):
    res = [1] * size
    for i in range(size):
        sum = 0
        for j in range(size):
            sum += mat[i][j] * vec[j]
        res[i] = sum
    return res

"""
    dot product
"""
def dot(vec1, vec2):
    size = len(vec1)
    s = 0
    for i in range(size):
        s += vec1[i] * vec2[i]
    return s

"""
    used after power series 
    for eigen vector to 
    normalize the eigen value
    (although scaing wont have any effect
    )
"""
def rayleigh_quotient(mat, vec):
    size = len(mat)
    matvec = multiply(mat, vec, size)
    return dot(matvec, vec) / dot(vec,vec)

def matrix_equal(mat1, mat2):
    size = len(mat1)
    equal = True
    precision = 0.1
    for i in mat1:
        for j in mat2:
            if abs( round(i -j, 3) ) > precision:
                return False
    return True

"""
    project the points onto the line 
    represented by this vector vec
"""
def project(l, vec):
    size = len(l)
    res = []
    for point in l:
        if point == (0,0):
            res.append( (0,0) )
            continue
        p = list(point)
        dot_product = dot(vec, p)
        mag1 = magnitude(vec)
        mag2 = magnitude(p)
        cosine = dot_product / (mag1 * mag2)
        proj_len = cosine * mag2

        proj_point = ( proj_len * vec[0] / mag1, proj_len * vec[1] / mag1)
        res.append(proj_point)
    return res

"""
    use power method
    gives principal eigen vector and value
"""
def eigen_problem(covar_mat):
    size = len(covar_mat)
    eigen_vector = [1] * size
    eigen_lambda = 1
    iterations = 100
    while True:
        iterations -= 1
        res, eigen_lambda = multiply_eigen(covar_mat, eigen_vector, size)
        eigen_vector = res
        if matrix_equal(eigen_vector, res) or iterations<1:
            break
    eigen_lambda = rayleigh_quotient(covar_mat, eigen_vector)
    return eigen_vector, eigen_lambda

"""
    return magnitude
"""
def magnitude(vector):
    s = 0
    for x in vector:
        s += x*x
    return round(math.sqrt(s), 3)

"""
    normalize the vector by magnitude
    unit vector
"""
def normalize_vector(vec):
    res = []
    size = len(vec)
    mag = magnitude(vec)
    for i in range(size):
        res.append( vec[i] / mag)
    return res

def point_plot(l, filename):
    x,y = zip(*l)
    plot.scatter(x,y)
    plot.show()
    plot.savefig(filename)
    plot.close()

def main():
    size = 3

    # randomly filled n points
    l = []
    mn = fill(l, 5, 10, 50)

    #l = [ (1,0), (2,0), (3,0), (4,4), (8,1) ]
    #l = [ (0,1), (0,3), (0,8), (4,4), (8,1) ]
    #mn = mean_calc(l)

    # original dataset
    point_plot(l, "original.png")

    normalize(l, mn)

    point_plot(l, "normalized.png")
    
    # covariance matrix
    covar_mat = covariance_matrix(l, mn)

    vector, value = eigen_problem(covar_mat)
    print(vector, value)
    vector = normalize_vector(vector)
    projected = project(l, vector)

    point_plot(projected, "projected.png")


if __name__=="__main__":
    main()

