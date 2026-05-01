module {
    func.func @load_guard(%arg0: i1, %arg1: i32, %arg2: i32) -> i32 {
        %res = rip.inline_mlir(%cond: i1 = %arg0, %ld: i32 = %arg1, %alt: i32 = %arg2) {
            %sel = arith.select %cond, %ld, %alt : i32
            rip.yield %sel : i32
        } -> i32
        return %res : i32
    }
}