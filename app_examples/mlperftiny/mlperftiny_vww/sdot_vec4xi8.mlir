module {
    func.func @sdot_vec4xi8(%arg0: i32, %arg1: i32) -> i32 {
        %res = rip.inline_mlir(%a: i32 = %arg0, %b: i32 = %arg1) {
            %a_vec = llvm.bitcast %a : i32 to vector<4xi8>
            %b_vec = llvm.bitcast %b : i32 to vector<4xi8>

            // On M+, this takes a total of 3 FUs, on Monaco it takes 4 FUs
            // M+: 1 ARITH/ARITH_PLUS, 2 BASIC_INT
            // Monaco: 2 FULL_INT, 1 BASIC_PLUS, 1 BASIC_INT

            // 1 ARITH/ARITH_PLUS/FULL_INT
            %low, %high = arith.mulsi_extended %a_vec, %b_vec: vector<4xi8>

            // 1 ARITH/ARITH_PLUS/FULL_INT - on M+ this gets merged into previous operation
            %0 = vector.shuffle %low, %high[0, 4, 1, 5]: vector<4xi8>, vector<4xi8>
            %1 = vector.shuffle %low, %high[2, 6, 3, 7]: vector<4xi8>, vector<4xi8>
            %cast0 = llvm.bitcast %0: vector<4xi8> to vector<2xi16>
            %cast1 = llvm.bitcast %1: vector<4xi8> to vector<2xi16>

            // Two additions - this is fast math
            // Note: -128*-128*2 = 32768 which does NOT fit in i16 but this
            // is the only case where overflow can happens so we accept it.
            // BASIC_INT
            %tmp = arith.addi %cast0, %cast1 : vector<2xi16>

            // Get's vector reduce
            // ARITH/ARITH_PLUS/BASIC_PLUS
            %tmp_low = vector.extract %tmp[0] : i16 from vector<2xi16>
            %tmp_high = vector.extract %tmp[1] : i16 from vector<2xi16>
            %tmp_low_signed = arith.extsi %tmp_low : i16 to i32
            %tmp_high_signed = arith.extsi %tmp_high : i16 to i32
            %sum = arith.addi %tmp_low_signed, %tmp_high_signed : i32

            rip.yield %sum : i32
        } -> i32
        return %res : i32
    }
}