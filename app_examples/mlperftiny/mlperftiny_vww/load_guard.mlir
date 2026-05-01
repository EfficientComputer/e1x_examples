module {
    func.func @load_guard(%arg0: i1, %arg1: i32, %arg2: i32) -> i32 {
        %res = rip.inline_mlir(%cond: i1 = %arg0, %ld: i32 = %arg1, %alt: i32 = %arg2) {
            %sel = arith.select %cond, %ld, %alt : i32
            rip.yield %sel : i32
        } -> i32
        return %res : i32
    }

    func.func @load_guard_i8(%arg0: i1, %arg1: i8, %arg2: i8) -> i8 {
        %res = rip.inline_mlir(
            %cond: i1 = %arg0,
            %ld:   i8 = %arg1,
            %alt:  i8 = %arg2
            ) {

            %sel = arith.select %cond, %ld, %alt : i8
            rip.yield %sel : i8
         } -> i8

        return %res : i8
    }

}