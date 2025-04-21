#
# M3L - Magayaga Mathematical Library (v0.9.2 / July 23, 2024)
# Copyright (c) 2024 Cyril John Magayaga (cjmagayaga957@gmail.com, cyrilmagayaga@proton.me)
#

class M3L
    PI = 3.141592653589793
    E = 2.718281828459045

    def self.positive(x)
        x
    end

    def self.negative(x)
        -x
    end
    
    def self.add(*args)
        result = 0
        args.each do |num|
            result += num
        end
        result
    end
  
    def self.subtract(*args)
        result = args[0]
        args[1..-1].each do |num|
            result -= num
        end
        result
    end
  
    def self.multiply(*args)
        result = 1
        args.each do |num|
            result *= num
        end
        result
    end
  
    def self.divide(*args)
        if args[1..-1].include?(0)
            raise ArgumentError, "Cannot divide by zero"
        end
        result = args[0]
        args[1..-1].each do |num|
            result /= num.to_f
        end
        result
    end

    def self.square(x)
        x ** 2
    end
    
    def self.cube(x)
        x ** 3
    end
    
    def self.sqrt(x)
        x ** (1/2)
    end
    
    def self.cbrt(x)
        x ** (1/3)
    end

    class Trim
        def self.sin(x)
            # Taylor series approximation for sin(x)
            result = 0
            10.times do |i|
                coef = (-1) ** i
                num = x ** (2 * i + 1)
                denom = M3L.factorial(2 * i + 1)
                result += coef * (num / denom)
            end
            result
        end
      
        def self.cos(x)
            # Taylor series approximation for cos(x)
            result = 0
            10.times do |i|
                coef = (-1) ** i
                num = x ** (2 * i)
                denom = self.factorial(2 * i)
                result += coef * (num / denom)
            end
            result
        end
      
        def self.tan(x)
            sin(x) / cos(x)
        end
      
        def self.csc(x)
            1 / sin(x)
        end
      
        def self.sec(x)
            1 / cos(x)
        end
      
        def self.cot(x)
            1 / tan(x)
        end

        def self.arcsin(x, terms = 10)
            raise ValueError, "Input should be in the range [-1, 1]" if x < -1 || x > 1
            
            result = x
            term = x
            x_squared = x * x
            (1...terms).each do |n|
              term *= x_squared * (2 * n - 1) / (2 * n).to_f
              result += term / (2 * n + 1).to_f
            end
            result
        end

        def self.arccos(x)
            raise ValueError, "Input should be in the range [-1, 1]" if x < -1 || x > 1
            
            # Abramowitz and Stegun approximation for arccos(x)
            M3L::PI / 2 - arcsin(x)
        end

        def self.arctan(z, terms = 10000)
            return 0 if z.zero?
        
            step = z.to_f / terms
            total_area = 0
            (0...terms).each do |i|
                t1 = i * step
                t2 = (i + 1) * step
                area = (1 / (1 + t1**2) + 1 / (1 + t2**2)) * step / 2.0
                total_area += area
            end
            total_area
        end
    end      

    def self.factorial(n)
        if n == 0
            1
        else
            n * factorial(n - 1)
        end
    end

    def self.log(base, x)
        raise ArgumentError, "Invalid input for logarithm" if base <= 0 || base == 1 || x <= 0
        ln(x) / ln(base)
    end

    def self.power(base, exp)
        base ** exp
    end

    def self.ln(x)
        raise ArgumentError, "Invalid input for natural logarithm" if x <= 0
        integrate(->(t) { 1 / t }, 1, x)
    end

    def self.integrate(f, a, b, n = 1000)
        dx = (b - a) / n
        integral = 0
        (0...n).each do |i|
            integral += f.call(a + (i + 0.5) * dx) * dx
        end
        integral
    end

    class Calc
        def self.summation(start, finish, term)
            raise ArgumentError, "Start index must be less than or equal to the end index" if start > finish
            result = 0
            (start..finish).each do |i|
                result += term.call(i)
            end
            result
        end
        
        def self.product(start, finish, term)
            raise ArgumentError, "Start index must be less than or equal to the end index" if start > finish
            result = 1
            (start..finish).each do |i|
                result *= term.call(i)
            end
            result
        end
    end

    def self.limit(func, x, approach="right", epsilon=1e-6)
        unless ["right", "left"].include?(approach)
            raise ValueError, "Approach must be 'right' or 'left'"
        end
        x_approach = approach == "right" ? x + epsilon : x - epsilon
        func.call(x_approach)
    end

    def self.array(matrix)
        matrix
    end

    def self.linspace(start, stop, num=50)
        raise ValueError, "Number of samples must be positive" if num <= 0
        step = (stop - start) / (num - 1)
        (0...num).map { |i| start + i * step }
    end

    def self.transpose(matrix)
        matrix.transpose
    end

    def self.inv(matrix)
        n = matrix.size
        raise ValueError, "Matrix must be square" if n != matrix[0].size
    
        identity = Array.new(n) { |i| Array.new(n, 0) }
        n.times { |i| identity[i][i] = 1 }
        
        n.times do |i|
            factor = 1.0 / matrix[i][i]
            n.times do |j|
                matrix[i][j] *= factor
                identity[i][j] *= factor
            end
            
            n.times do |k|
                next if k == i
                factor = matrix[k][i]
                n.times do |j|
                    matrix[k][j] -= factor * matrix[i][j]
                    identity[k][j] -= factor * identity[i][j]
                end
            end
        end
        identity
    end

    def self.dot(a, b)
        if a[0].is_a?(Array) && b[0].is_a?(Array)
            a.size.times.map do |i|
                b.transpose.map { |col| col.each_with_index.map { |elem, j| elem * a[i][j] }.sum }
            end
        else
            a.zip(b).map { |x, y| x * y }.sum
        end
    end

    def self.vdot(a, b)
        if a[0].is_a?(Array) && b[0].is_a?(Array)
            a.flatten.zip(b.flatten).map { |x, y| x * y }.sum
        else
            a.zip(b).map { |x, y| x * y }.sum
        end
    end

    def self.multi_dot(matrices)
        raise ValueError, "At least two matrices are required for multi_dot" if matrices.size < 2
    
        result = matrices[0]
        matrices[1..-1].each { |matrix| result = multiply_matrices(result, matrix) }
        result
    end

    def self.random(size)
        if size.is_a?(Integer)
            Array.new(size) { rand }
        elsif size.is_a?(Array)
            size.size == 1 ? Array.new(size[0]) { rand } : Array.new(size[0]) { Array.new(size[1]) { rand } }
        else
            raise ValueError, "Invalid size for random"
        end
    end

    def self.solve(matrix, b)
        inv_matrix = inv(matrix)
        inv_matrix.map { |row| row.zip(b).map { |m, n| m * n }.sum }
    end

    def self.multiply_matrices(matrix1, matrix2)
        raise ValueError, "Matrix dimensions incompatible for multiplication" if matrix1[0].size != matrix2.size
    
        result = Array.new(matrix1.size) { Array.new(matrix2[0].size, 0) }
        matrix1.size.times do |i|
            matrix2[0].size.times do |j|
                matrix1[0].size.times do |k|
                    result[i][j] += matrix1[i][k] * matrix2[k][j]
                end
            end
        end
        result
    end

    def self.matrix_power(matrix, n)
        return Array.new(matrix.size) { |i| Array.new(matrix.size) { |j| (i == j) ? 1 : 0 } } if n == 0
        return matrix if n == 1
    
        result = matrix
        (n - 1).times { result = multiply_matrices(result, matrix) }
        result
    end

    def self.arange(start, stop = nil, step = 1)
        stop.nil? ? (0...start).step(step).to_a : (start...stop).step(step).to_a
    end

    def self.reshape(matrix, new_shape)
        flat = matrix.flatten
        reshaped = []
        count = 0
        new_shape[0].times do |i|
            row = []
            new_shape[1].times do |j|
                row << flat[count]
                count += 1
            end
            reshaped << row
        end
        reshaped
    end

    def self.zeros(shape)
        case shape
        when Integer then Array.new(shape, 0)
        when Array
            shape.size == 1 ? Array.new(shape[0], 0) : Array.new(shape[0]) { Array.new(shape[1], 0) }
        else raise ValueError, "Invalid shape for zeros"
        end
    end

    def self.tensordot(a, b, axes = 2)
        if a[0].is_a?(Array) && b[0].is_a?(Array)
            a.flatten.zip(b.flatten).map { |x, y| x * y }.sum
        else
            a.zip(b).map { |x, y| x * y }.sum
        end
    end

    def self.kronecker(matrix1, matrix2)
        result = []
        matrix1.each do |row1|
            matrix2.each do |row2|
                result_row = row1.map { |x| x * row2 }
                result << result_row
            end
        end
        result
    end

    def self.cholesky(matrix)
        n = matrix.size
        l = Array.new(n) { Array.new(n, 0) }
    
        n.times do |i|
            (i + 1).times do |j|
                sum = 0
                if i == j
                    (0...j).each { |k| sum += l[j][k] ** 2 }
                    l[j][j] = Math.sqrt(matrix[j][j] - sum)
                else
                    (0...j).each { |k| sum += l[i][k] * l[j][k] }
                    l[i][j] = (matrix[i][j] - sum) / l[j][j]
                end
            end
        end
        l
    end

    class Binary
        def self.bdn(binary)
            binary.to_s.chars.reverse.each_with_index.map { |bit, i| bit.to_i * (2 ** i) }.sum
        end
    
        def self.dbn(decimal)
            return "0" if decimal.zero?
        
            binary = ""
            while decimal > 0
                binary.prepend((decimal % 2).to_s)
                decimal /= 2
            end
            binary
        end
    end
    
    def self.prime_factorization(number)
        factors = []
        divisor = 2
        while number > 1
            if number % divisor == 0
                factors << divisor
                number /= divisor
            else
                divisor += 1
            end
        end
        factors.join(" * ")
    end

    class Function
        def self.gamma(z)
            sqrt_two_pi = self.sqrt(M2L::PI * 2)
            sqrt_two_pi * ((z - 1/2) ** (z - 1/2)) / (self::E ** z)
        end
    
        def self.zeta(s)
            result = 0
            n = 1
            loop do
                term = 1 / (n ** s)
                break if term < 1e-6
                result += term
                n += 1
            end
            result
        end
    end

end