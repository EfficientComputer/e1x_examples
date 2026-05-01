module {
    func.func @add_tree_8(
        %arg0: i32, %arg1: i32, %arg2: i32, %arg3: i32,
        %arg4: i32, %arg5: i32, %arg6: i32, %arg7: i32
    ) -> i32 {
        %res = rip.inline_mlir(
            %a0: i32 = %arg0, %a1: i32 = %arg1,
            %a2: i32 = %arg2, %a3: i32 = %arg3,
            %a4: i32 = %arg4, %a5: i32 = %arg5,
            %a6: i32 = %arg6, %a7: i32 = %arg7
        ) {

            // Build an addition tree
            %sum0 = arith.addi %a0, %a1 : i32
            %sum1 = arith.addi %a2, %a3 : i32
            %sum2 = arith.addi %a4, %a5 : i32
            %sum3 = arith.addi %a6, %a7 : i32

            %sum4 = arith.addi %sum0, %sum1 : i32
            %sum5 = arith.addi %sum2, %sum3 : i32

            %total = arith.addi %sum4, %sum5 : i32

            rip.yield %total : i32
        } -> i32
        return %res : i32
    }

    func.func @add_tree_10(
        %arg0: i32, %arg1: i32, %arg2: i32, %arg3: i32, %arg4: i32,
        %arg5: i32, %arg6: i32, %arg7: i32, %arg8: i32, %arg9: i32
    ) -> i32 {
        %res = rip.inline_mlir(
            %a0: i32 = %arg0, %a1: i32 = %arg1, %a2: i32 = %arg2, %a3: i32 = %arg3, %a4: i32 = %arg4,
            %a5: i32 = %arg5, %a6: i32 = %arg6, %a7: i32 = %arg7, %a8: i32 = %arg8, %a9: i32 = %arg9
        ) {
            // Pairwise
            %s0 = arith.addi %a0, %a1 : i32
            %s1 = arith.addi %a2, %a3 : i32
            %s2 = arith.addi %a4, %a5 : i32
            %s3 = arith.addi %a6, %a7 : i32
            %s4 = arith.addi %a8, %a9 : i32

            // Reduce (balanced-ish)
            %t0 = arith.addi %s0, %s1 : i32
            %t1 = arith.addi %s2, %s3 : i32

            %u0 = arith.addi %t0, %t1 : i32
            %total = arith.addi %u0, %s4 : i32

            rip.yield %total : i32
        } -> i32
        return %res : i32
    }
}