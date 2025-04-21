--
-- M3L - Magayaga Mathematical Library (v0.9.2 / July 23, 2024)
-- Copyright (c) 2024-2025 Cyril John Magayaga (cjmagayaga957@gmail.com, cyrilmagayaga@proton.me)
--

-- Define the M3L class
M3L = {}

-- M3L Trigonometry
M3L.Trim = {}

-- M3L Calculus
M3L.Calc = {}

-- M3L Binary
M3L.Binary = {}

-- M3L Math Function
M3L.Function = {}

-- Define class variables
M3L.PI = 3.141592653589793
M3L.TAU = 2 * M3L.PI
M3L.E = 2.718281828459045

-- Define static methods
function M3L.positive(x)
    return x
end

function M3L.negative(x)
    return -x
end

function M3L.add(...)
    local args = {...}
    local result = 0
    for _, num in ipairs(args) do
        result = result + num
    end
    return result
end

function M3L.subtract(...)
    local args = {...}
    local result = args[1]
    for i = 2, #args do
        result = result - args[i]
    end
    return result
end

function M3L.multiply(...)
    local args = {...}
    local result = 1
    for _, num in ipairs(args) do
        result = result * num
    end
    return result
end

function M3L.divide(...)
    local args = {...}
    if #args < 2 then
        error("At least two arguments are required for division")
    end
    for _, num in ipairs(args) do
        if num == 0 then
            error("Cannot divide by zero")
        end
    end
    local result = args[1]
    for i = 2, #args do
        result = result / args[i]
    end
    return result
end

function M3L.square(x)
    return x ^ 2
end

function M3L.cube(x)
    return x ^ 3
end

function M3L.sqrt(x)
    return x ^ (1/2)
end

function M3L.cbrt(x)
    return x ^ (1/3)
end

function M3L.abs(x)
    return x < 0 and -x or x
end

function M3L.gcd(a, b)
    while b ~= 0 do
        a, b = b, a % b
    end
    return M3L.abs(a)
end

function M3L.lcm(a, b)
    if a == 0 or b == 0 then
        return 0
    end
    return M3L.abs(a * b) // M3L.gcd(a, b)
end

function M3L.fma(a, b, c)
    return a * b + c
end

function M3L.floor(x)
    local i = x // 1  -- Truncate towards zero
    if x < 0 and x ~= i then
        return i - 1
    else
        return i
    end
end

function M3L.Trim.sin(x)
    -- Taylor series approximation for sin(x)
    local result = 0
    for i = 0, 9 do
        local coef = (-1) ^ i
        local num = x ^ (2 * i + 1)
        local denom = M3L.factorial(2 * i + 1)
        result = result + coef * (num / denom)
    end
    return result
end

function M3L.Trim.cos(x)
    -- Taylor series approximation for cos(x)
    local result = 0
    for i = 0, 9 do
        local coef = (-1) ^ i
        local num = x ^ (2 * i)
        local denom = M3L.factorial(2 * i)
        result = result + coef * (num / denom)
    end
    return result
end

function M3L.Trim.tan(x)
    return M3L.sin(x) / M3L.cos(x)
end

function M3L.Trim.csc(x)
    return 1 / M3L.sin(x)
end

function M3L.Trim.sec(x)
    return 1 / M3L.cos(x)
end

function M3L.Trim.cot(x)
    return 1 / M3L.tan(x)
end

function M3L.Trim.arcsin(x, terms)
    terms = terms or 10
    if x < -1 or x > 1 then
        error("Input should be in the range [-1, 1]", 2)
    end
    result = x
    term = x
    x_squared = x * x
    for n = 1, terms - 1 do
        term = term * x_squared * (2 * n - 1) / (2 * n)
        result = result + term / (2 * n + 1)
    end
    return result
end

function M3L.Trim.arccos(x)
    if x < -1 or x > 1 then
        error("Input should be in the range [-1, 1]", 2)
    end
    -- Abramowitz and Stegun approximation for arccos(x)
    return M3L.PI / 2 - M3L.Trim.arcsin(x)
end

function M3L.Trim.arctan(z, terms)
    terms = terms or 10000
    if z == 0 then
        return 0
    end
    step = z / terms
    total_area = 0
    for i = 0, terms - 1 do
        t1 = i * step
        t2 = (i + 1) * step
        area = (1 / (1 + t1^2) + 1 / (1 + t2^2)) * step / 2
        total_area = total_area + area
    end
    return total_area
end

function M3L.factorial(n)
    if n == 0 then
        return 1
    else
        return n * M3L.factorial(n - 1)
    end
end

function M3L.log(base, x)
    if base <= 0 or base == 1 or x <= 0 then
        error("Invalid input for logarithm")
    end
    return M3L.ln(x) / M3L.ln(base)
end

function M3L.power(base, exp)
    return base ^ exp
end

function M3L.ln(x)
    if x <= 0 then
        error("Invalid input for natural logarithm")
    end
    return M3L.integrate(function(t) return 1 / t end, 1, x)
end

function M3L.Calc.summation(start, finish, term)
    if start > finish then
        error("Start index must be less than or equal to the end index")
    end
    local result = 0
    for i = start, finish do
        result = result + term(i)
    end
    return result
end

function M3L.Calc.product(start, finish, term)
    if start > finish then
        error("Start index must be less than or equal to the end index")
    end
    local result = 1
    for i = start, finish do
        result = result * term(i)
    end
    return result
end

function M3L.integrate(f, a, b, N)
    N = N or 1000
    local dx = (b - a) / N
    local integral = 0
    for i = 0, N - 1 do
        integral = integral + f(a + (i + 0.5) * dx) * dx
    end
    return integral
end

function M3L.limit(func, x, approach, epsilon)
    approach = approach or "right"
    epsilon = epsilon or 1e-6
    if approach ~= "right" and approach ~= "left" then
        error("Approach must be 'right' or 'left'")
    end
    local x_approach
    if approach == "right" then
        x_approach = x + epsilon
    else
        x_approach = x - epsilon
    end
    return func(x_approach)
end

function M3L.array(matrix)
    return matrix
end

function M3L.linspace(start, stop, num)
    num = num or 50
    if num <= 0 then
        error("Number of samples must be positive")
    end
    local step = (stop - start) / (num - 1)
    local result = {}
    for i = 0, num - 1 do
        table.insert(result, start + i * step)
    end
    return result
end

function M3L.transpose(matrix)
    local result = {}
    for i = 1, #matrix[1] do
        local row = {}
        for j = 1, #matrix do
            table.insert(row, matrix[j][i])
        end
        table.insert(result, row)
    end
    return result
end

function M3L.inv(matrix)
    local n = #matrix
    if n ~= #matrix[1] then
        error("Matrix must be square")
    end
    local identity = {}
    for i = 1, n do
        identity[i] = {}
        for j = 1, n do
            identity[i][j] = 0
        end
        identity[i][i] = 1
    end
    for i = 1, n do
        local factor = 1 / matrix[i][i]
        for j = 1, n do
            matrix[i][j] = matrix[i][j] * factor
            identity[i][j] = identity[i][j] * factor
        end
        for k = 1, n do
            if k ~= i then
                factor = matrix[k][i]
                for j = 1, n do
                    matrix[k][j] = matrix[k][j] - factor * matrix[i][j]
                    identity[k][j] = identity[k][j] - factor * identity[i][j]
                end
            end
        end
    end
    return identity
end

function M3L.dot(a, b)
    if type(a[1]) == "table" and type(b[1]) == "table" then
        local result = {}
        for i = 1, #a do
            result[i] = {}
            for j = 1, #b[1] do
                local sum = 0
                for k = 1, #a[1] do
                    sum = sum + a[i][k] * b[k][j]
                end
                result[i][j] = sum
            end
        end
        return result
    else
        local sum = 0
        for i = 1, #a do
            sum = sum + a[i] * b[i]
        end
        return sum
    end
end

function M3L.vdot(a, b)
    if type(a[1]) == "table" and type(b[1]) == "table" then
        local a_flat = {}
        for _, sublist in ipairs(a) do
            for _, elem in ipairs(sublist) do
                table.insert(a_flat, elem)
            end
        end
        local b_flat = {}
        for _, sublist in ipairs(b) do
            for _, elem in ipairs(sublist) do
                table.insert(b_flat, elem)
            end
        end
        local sum = 0
        for i = 1, #a_flat do
            sum = sum + a_flat[i] * b_flat[i]
        end
        return sum
    else
        local sum = 0
        for i = 1, #a do
            sum = sum + a[i] * b[i]
        end
        return sum
    end
end

function M3L.multi_dot(matrices)
    if #matrices < 2 then
        error("At least two matrices are required for multi_dot")
    end
    local result = matrices[1]
    for i = 2, #matrices do
        result = M3L.dot(result, matrices[i])
    end
    return result
end

function M3L.random(size)
    if type(size) == "number" then
        local result = {}
        for _ = 1, size do
            table.insert(result, math.random())
        end
        return result
    elseif type(size) == "table" then
        if #size == 1 then
            local result = {}
            for _ = 1, size[1] do
                table.insert(result, math.random())
            end
            return result
        elseif #size == 2 then
            local result = {}
            for _ = 1, size[1] do
                local row = {}
                for _ = 1, size[2] do
                    table.insert(row, math.random())
                end
                table.insert(result, row)
            end
            return result
        else
            error("Invalid size for random")
        end
    else
        error("Invalid size for random")
    end
end


function M3L.solve(matrix, b)
    local inv_matrix = M3L.inv(matrix)
    local result = {}
    for i = 1, #matrix do
        local sum = 0
        for j = 1, #b do
            sum = sum + inv_matrix[i][j] * b[j]
        end
        table.insert(result, sum)
    end
    return result
end

function M3L.multiply_matrices(matrix1, matrix2)
    local result = {}
    for i = 1, #matrix1 do
        result[i] = {}
        for j = 1, #matrix2[1] do
            local sum = 0
            for k = 1, #matrix2 do
                sum = sum + matrix1[i][k] * matrix2[k][j]
            end
            result[i][j] = sum
        end
    end
    return result
end

function M3L.matrix_power(matrix, n)
    if n == 0 then
        local result = {}
        for i = 1, #matrix do
            result[i] = {}
            for j = 1, #matrix do
                result[i][j] = i == j and 1 or 0
            end
        end
        return result
    elseif n > 0 then
        local result = matrix
        for _ = 2, n do
            result = M3L.multiply_matrices(result, matrix)
        end
        return result
    else
        error("Matrix power is only defined for non-negative integers")
    end
end

function M3L.arange(start, stop, step)
    if stop == nil then
        start, stop = 0, start
    end
    step = step or 1
    local result = {}
    for i = start, stop - 1, step do
        table.insert(result, i)
    end
    return result
end

function M3L.reshape(matrix, new_shape)
    local flat = {}
    for _, sublist in ipairs(matrix) do
        for _, item in ipairs(sublist) do
            table.insert(flat, item)
        end
    end
    local reshaped = {}
    local count = 1
    for i = 1, new_shape[1] do
        local row = {}
        for j = 1, new_shape[2] do
            row[j] = flat[count]
            count = count + 1
        end
        reshaped[i] = row
    end
    return reshaped
end

function M3L.zeros(shape)
    if type(shape) == "number" then
        return {0}
    elseif type(shape) == "table" then
        if #shape == 1 then
            local result = {}
            for _ = 1, shape[1] do
                table.insert(result, 0)
            end
            return result
        elseif #shape == 2 then
            local result = {}
            for _ = 1, shape[1] do
                local row = {}
                for _ = 1, shape[2] do
                    table.insert(row, 0)
                end
                table.insert(result, row)
            end
            return result
        else
            error("Invalid shape for zeros")
        end
    else
        error("Invalid shape for zeros")
    end
end

function M3L.tensordot(a, b, axes)
    if type(a[1]) == "table" and type(b[1]) == "table" then
        local a_flat = {}
        for _, sublist in ipairs(a) do
            for _, item in ipairs(sublist) do
                table.insert(a_flat, item)
            end
        end
        local b_flat = {}
        for _, sublist in ipairs(b) do
            for _, item in ipairs(sublist) do
                table.insert(b_flat, item)
            end
        end
        local sum = 0
        for i = 1, #a_flat do
            sum = sum + a_flat[i] * b_flat[i]
        end
        return sum
    else
        local sum = 0
        for i = 1, #a do
            sum = sum + a[i] * b[i]
        end
        return sum
    end
end

function M3L.kronecker(matrix1, matrix2)
    local result = {}
    for i = 1, #matrix1 do
        for j = 1, #matrix1[1] do
            for m = 1, #matrix2 do
                for n = 1, #matrix2[1] do
                    if not result[i] then result[i] = {} end
                    if not result[i][j] then result[i][j] = {} end
                    table.insert(result[i][j], matrix1[i][j] * matrix2[m][n])
                end
            end
        end
    end
    return result
end

function M3L.cholesky(matrix)
    local n = #matrix
    local L = {}
    for i = 1, n do
        L[i] = {}
        for j = 1, n do
            L[i][j] = 0
        end
    end
    for i = 1, n do
        for j = 1, i do
            if i == j then
                local sum_sq = 0
                for k = 1, j - 1 do
                    sum_sq = sum_sq + L[i][k] ^ 2
                end
                L[i][j] = math.sqrt(matrix[i][i] - sum_sq)
            else
                local sum_prod = 0
                for k = 1, j - 1 do
                    sum_prod = sum_prod + L[i][k] * L[j][k]
                end
                L[i][j] = (matrix[i][j] - sum_prod) / L[j][j]
            end
        end
    end
    return L
end

function M3L.Binary.bdn(binary)
    local binary_str = tostring(binary)
    local decimal = 0
    for i = 1, #binary_str do
        local bit = tonumber(binary_str:sub(i, i))
        if bit ~= 0 and bit ~= 1 then
            error("Invalid binary number")
        end
        decimal = decimal + bit * (2 ^ (#binary_str - i))
    end
    return decimal
end

function M3L.Binary.dbn(decimal)
    if decimal < 0 then
        error("Decimal number must be non-negative")
    end
    if decimal == 0 then
        return "0"
    end
    local binary = ""
    while decimal > 0 do
        binary = tostring(decimal % 2) .. binary
        decimal = math.floor(decimal / 2)
    end
    return binary
end

function M3L.prime_factorization(number)
    local factors = {}
    local divisor = 2
    while number > 1 do
        if number % divisor == 0 then
            table.insert(factors, divisor)
            number = number / divisor
        else
            divisor = divisor + 1
        end
    end
    return table.concat(factors, " * ")
end

function M3L.Function.gamma(z)
    local sqrt_two_pi = math.sqrt(M3L.PI * 2)
    return sqrt_two_pi * ((z - 1/2) ^ (z - 1/2)) / (M3L.E ^ z)
end

function M3L.Function.zeta(s)
    local result = 0
    local n = 1
    while true do
        local term = 1 / (n ^ s)
        if term < 1e-6 then
            break
        end
        result = result + term
        n = n + 1
    end
    return result
end