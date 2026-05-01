module {
    func.func @load_guard_i16(%arg0: i1, %arg1: i16, %arg2: i16) -> i16 {
        %res = rip.inline_mlir(%cond: i1 = %arg0, %ld: i16 = %arg1, %alt: i16 = %arg2) {
            %sel = arith.select %cond, %ld, %alt : i16
            rip.yield %sel : i16
        } -> i16
        return %res : i16
    }
    func.func @load_guard_i32(%arg0: i1, %arg1: i32, %arg2: i32) -> i32 {
        %res = rip.inline_mlir(%cond: i1 = %arg0, %ld: i32 = %arg1, %alt: i32 = %arg2) {
            %sel = arith.select %cond, %ld, %alt : i32
            rip.yield %sel : i32
        } -> i32
        return %res : i32
    }
}