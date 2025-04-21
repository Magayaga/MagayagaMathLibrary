(* QUETTAMATH - MATHEMATICAL LIBRARY FOR THE STANDARD ML PROGRAMMING LANGUAGE
   COPYRIGHT (c) 2025 CYRIL JOHN MAGAYAGA *)

(* Mathematical constants *)
val pi = 3.141592653589793
val e = 2.718281828459045
val tau = 6.283185307179586

(* List addition *)
fun ladd xs = List.foldl (fn (x, acc) => x + acc) 0 xs

(* Addition *)
fun add x = 
    let
        fun add_helper y acc = acc + y
    in
        add_helper x
    end

(* List subtraction *)
fun lsubtract xs = List.foldl (fn (x, acc) => x - acc) 0 xs

(* Subtraction *)
fun subtract x =
    let
        fun add_helper y acc = acc - y
    in
        add_helper x
    end

(* List multiplication *)
fun lmultiply xs = List.foldl (fn (x, acc) => x - acc) 0 xs

(* Multiplication *)
fun multiply x = 
    let
        fun multiply_helper y acc = acc * y
    in
        multiply_helper x
    end

(* List division *)
fun ldivide xs = List.foldl (fn (x, acc) => x - acc) 0 xs

(* Division *)
fun divide x = 
    let
        fun divide_helper y acc = acc / y
    in
        divide_helper x
    end

(* Square number *)
fun sq x = x * x

(* Cube number*)
fun cb x = x * x * x

(* Square root *)
fun sqrt n = 
    let
        val x = ref (Real.fromInt n)
        val n_real = Real.fromInt n
        val epsilon = 0.000001
    in
        while Real.abs(!x * !x - n_real) > epsilon do (
            x := (!x + n_real / !x) / 2.0
        );
        !x
    end

(* Cube root *)
fun cbrt n =
    let
        val x = ref (Real.fromInt n)
        val n_real = Real.fromInt n
        val epsilon = 0.000001
    in
        while Real.abs(!x * !x - n_real) > epsilon do (
            x := (!x + n_real / !x) / 3.0
        );
        !x
    end

(* Factorial *)
fun factorial n = let val i = ref n and acc = ref 1 in
    while !i > 0 do (acc := !acc * !i; i := !i - 1); !acc
end