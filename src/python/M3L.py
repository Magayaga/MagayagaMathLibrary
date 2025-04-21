#
# M3L - Magayaga Mathematical Library (v0.9.2 / July 23, 2024)
# Copyright (c) 2024 Cyril John Magayaga (cjmagayaga957@gmail.com, cyrilmagayaga@proton.me)
#

import random

class M3L:
    PI = 3.141592653589793
    E = 2.718281828459045

    @staticmethod
    def positive(x):
        return x

    @staticmethod
    def negative(x):
        return -x
    
    @staticmethod
    def add(*args):
        result = 0
        for num in args:
            result += num
        return result

    @staticmethod
    def subtract(*args):
        result = args[0]
        for num in args[1:]:
            result -= num
        return result

    @staticmethod
    def multiply(*args):
        result = 1
        for num in args:
            result *= num
        return result

    @staticmethod
    def divide(*args):
        if 0 in args[1:]:
            raise ValueError("Cannot divide by zero")
        result = args[0]
        for num in args[1:]:
            result /= num
        return result
    
    @staticmethod
    def square(x):
        return x ** 2

    @staticmethod
    def cube(x):
        return x ** 3

    @staticmethod
    def sqrt(x):
        return x ** (1/2)

    @staticmethod
    def cbrt(x):
        return x ** (1/3)
    
    class trim:
        @staticmethod
        def sin(x):
            # Taylor series approximation for sin(x)
            result = 0
            for i in range(10):
                coef = (-1) ** i
                num = x ** (2 * i + 1)
                denom = M3L.factorial(2 * i + 1)
                result += coef * (num / denom)
            return result
    
        @staticmethod
        def cos(x):
            # Taylor series approximation for cos(x)
            result = 0
            for i in range(10):
                coef = (-1) ** i
                num = x ** (2 * i)
                denom = M3L.factorial(2 * i)
                result += coef * (num / denom)
            return result
    
        @staticmethod
        def tan(x):
            return M3L.sin(x) / M3L.cos(x)
    
        @staticmethod
        def csc(x):
            return 1 / M3L.sin(x)
    
        @staticmethod
        def sec(x):
            return 1 / M3L.cos(x)
    
        @staticmethod
        def cot(x):
            return 1 / M3L.tan(x)
        
        @staticmethod
        def arcsin(x, terms=10):
            if x < -1 or x > 1:
                raise ValueError("Input should be in the range [-1, 1]")
            result = x
            term = x
            x_squared = x * x
            for n in range(1, terms):
                term *= x_squared * (2 * n - 1) / (2 * n)
                result += term / (2 * n + 1)
            return result

        @staticmethod
        def arccos(x):
            if x < -1 or x > 1:
                raise ValueError("Input should be in the range [-1, 1]")
            # Abramowitz and Stegun approximation for arccos(x)
            return M3L.PI / 2 - M3L.trim.arcsin(x)
        
        @staticmethod
        def arctan(z, terms=10000):
            if z == 0:
                return 0
            step = z / terms
            total_area = 0
            for i in range(terms):
                t1 = i * step
                t2 = (i + 1) * step
                area = (1 / (1 + t1**2) + 1 / (1 + t2**2)) * step / 2
                total_area += area
            return total_area

    @staticmethod
    def factorial(n):
        if n == 0:
            return 1
        else:
            return n * M3L.factorial(n - 1)
    
    @staticmethod
    def log(base, x):
        if base <= 0 or base == 1 or x <= 0:
            raise ValueError("Invalid input for logarithm")
        return M3L.ln(x) / M3L.ln(base)

    @staticmethod
    def power(base, exp):
        return base ** exp

    @staticmethod
    def ln(x):
        if x <= 0:
            raise ValueError("Invalid input for natural logarithm")
        return M3L.integrate(lambda t: 1 / t, 1, x)
    
    class calc:
        @staticmethod
        def summation(start, end, term):
            if start > end:
                raise ValueError("Start index must be less than or equal to the end index")
            result = 0
            for i in range(start, end + 1):
                result += term(i)
            return result
        
        @staticmethod
        def product(start, end, term):
            if start > end:
                raise ValueError("Start index must be less than or equal to the end index")
            result = 1
            for i in range(start, end + 1):
                result *= term(i)
            return result
        
    @staticmethod
    def integrate(f, a, b, N=1000):
        dx = (b - a) / N
        integral = 0
        for i in range(N):
            integral += f(a + (i + 0.5) * dx) * dx
        return integral
    
    # Example code:
    #
    # You can be f(x), f(x) = M3L.sin(x).
    #
    # result = M3L.limit(f, math.pi/2, approach='right')
    # print(result)

    @staticmethod
    def limit(func, x, approach="right", epsilon=1e-6):
        if approach not in ("right", "left"):
            raise ValueError("Approach must be 'right' or 'left'")
        if approach == "right":
            x_approach = x + epsilon
        else:
            x_approach = x - epsilon
        return func(x_approach)

    @staticmethod
    def array(matrix):
        return matrix

    @staticmethod
    def linspace(start, stop, num=50):
        if num <= 0:
            raise ValueError("Number of samples must be positive")
        step = (stop - start) / (num - 1)
        return [start + i * step for i in range(num)]

    @staticmethod
    def transpose(matrix):
        return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]
    
    @staticmethod
    def inv(matrix):
        n = len(matrix)
        if n != len(matrix[0]):
            raise ValueError("Matrix must be square")
        identity = [[0] * n for _ in range(n)]
        for i in range(n):
            identity[i][i] = 1
        for i in range(n):
            factor = 1 / matrix[i][i]
            for j in range(n):
                matrix[i][j] *= factor
                identity[i][j] *= factor
            for k in range(n):
                if k != i:
                    factor = matrix[k][i]
                    for j in range(n):
                        matrix[k][j] -= factor * matrix[i][j]
                        identity[k][j] -= factor * identity[i][j]
        return identity

    @staticmethod
    def dot(a, b):
        if isinstance(a[0], (list, tuple)) and isinstance(b[0], (list, tuple)):
            return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]
        else:
            return sum(x * y for x, y in zip(a, b))
    
    @staticmethod
    def vdot(a, b):
        if isinstance(a[0], (list, tuple)) and isinstance(b[0], (list, tuple)):
            a_flat = [elem for sublist in a for elem in sublist]
            b_flat = [elem for sublist in b for elem in sublist]
            return sum(x * y for x, y in zip(a_flat, b_flat))
        else:
            return sum(x * y for x, y in zip(a, b))

    @staticmethod
    def multi_dot(matrices):
        if len(matrices) < 2:
            raise ValueError("At least two matrices are required for multi_dot")
        result = matrices[0]
        for matrix in matrices[1:]:
            result = M3L.multiply_matrices(result, matrix)
        return result

    @staticmethod
    def random(size):
        if isinstance(size, int):
            return [random.random() for _ in range(size)]
        elif isinstance(size, tuple):
            if len(size) == 1:
                return [random.random() for _ in range(size[0])]
            elif len(size) == 2:
                return [[random.random() for _ in range(size[1])] for _ in range(size[0])]
            else:
                raise ValueError("Invalid size for random")
        else:
            raise ValueError("Invalid size for random")

    @staticmethod
    def solve(matrix, b):
        inv_matrix = M3L.inv(matrix)
        return [sum(inv_matrix[i][j] * b[j] for j in range(len(b))) for i in range(len(matrix))]
    
    @staticmethod
    def multiply_matrices(matrix1, matrix2):
        if len(matrix1[0]) != len(matrix2):
            raise ValueError("Matrix dimensions incompatible for multiplication")
        result = [[0] * len(matrix2[0]) for _ in range(len(matrix1))]
        for i in range(len(matrix1)):
            for j in range(len(matrix2[0])):
                for k in range(len(matrix2)):
                    result[i][j] += matrix1[i][k] * matrix2[k][j]
        return result
    
    @staticmethod
    def matrix_power(matrix, n):
        if n == 0:
            return [[1 if i == j else 0 for j in range(len(matrix))] for i in range(len(matrix))]
        elif n > 0:
            result = matrix
            for _ in range(n - 1):
                result = M3L.multiply_matrices(result, matrix)
            return result
        else:
            raise ValueError("Matrix power is only defined for non-negative integers")
    
    @staticmethod
    def arange(start, stop=None, step=1):
        if stop is None:
            start, stop = 0, start
        return [i for i in range(start, stop, step)]

    @staticmethod
    def reshape(matrix, new_shape):
        flat = [item for sublist in matrix for item in sublist]
        reshaped = []
        count = 0
        for i in range(new_shape[0]):
            row = []
            for j in range(new_shape[1]):
                row.append(flat[count])
                count += 1
            reshaped.append(row)
        return reshaped

    @staticmethod
    def zeros(shape):
        if isinstance(shape, int):
            return [0] * shape
        elif isinstance(shape, tuple):
            if len(shape) == 1:
                return [0] * shape[0]
            elif len(shape) == 2:
                return [[0] * shape[1] for _ in range(shape[0])]
            else:
                raise ValueError("Invalid shape for zeros")
        else:
            raise ValueError("Invalid shape for zeros")

    @staticmethod
    def tensordot(a, b, axes=2):
        if isinstance(a[0], list) and isinstance(b[0], list):
            a_flat = [elem for sublist in a for elem in sublist]
            b_flat = [elem for sublist in b for elem in sublist]
            return sum(x * y for x, y in zip(a_flat, b_flat))
        else:
            return sum(x * y for x, y in zip(a, b))

    @staticmethod
    def kronecker(matrix1, matrix2):
        result = []
        for i in range(len(matrix1)):
            row = []
            for j in range(len(matrix1[0])):
                for m in range(len(matrix2)):
                    for n in range(len(matrix2[0])):
                        row.append(matrix1[i][j] * matrix2[m][n])
            result.append(row)
        return result

    @staticmethod
    def cholesky(matrix):
        n = len(matrix)
        L = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1):
                if i == j:
                    sum_sq = sum(L[i][k] ** 2 for k in range(j))
                    L[i][j] = M3L.sqrt(matrix[i][i] - sum_sq)
                else:
                    sum_prod = sum(L[i][k] * L[j][k] for k in range(j))
                    L[i][j] = (matrix[i][j] - sum_prod) / L[j][j]
        return L
    
    class binary:
        @staticmethod
        def bdn(binary):
            decimal = 0
            binary_str = str(binary)
            for i in range(len(binary_str)):
                bit = int(binary_str[i])
                if bit != 0 and bit != 1:
                    raise ValueError("Invalid binary number")
                decimal += bit * (2 ** (len(binary_str) - i - 1))
            return decimal
        
        @staticmethod
        def dbn(decimal):
            if decimal < 0:
                raise ValueError("Decimal number must be non-negative")
            binary = ""
            if decimal == 0:
                return "0"
            while decimal > 0:
                binary = str(decimal % 2) + binary
                decimal //= 2
            return binary
        
    @staticmethod
    def prime_factorization(number):
        factors = []
        divisor = 2
        while number > 1:
            if number % divisor == 0:
                factors.append(divisor)
                number //= divisor
            else:
                divisor += 1
        return " * ".join(str(factor) for factor in factors)
    
    class function:
        @staticmethod
        def gamma(z):
            sqrt_two_pi = M3L.sqrt(M3L.PI * 2)
            return sqrt_two_pi * ((z - 1/2) ** (z - 1/2)) / (M3L.E ** z)
        
        @staticmethod 
        def zeta(s):
            result = 0
            n = 1
            while True:
                term = 1 / (n ** s)
                if term < 1e-6:  # Adjust this threshold for desired precision
                    break
                result += term
                n += 1
                return result
