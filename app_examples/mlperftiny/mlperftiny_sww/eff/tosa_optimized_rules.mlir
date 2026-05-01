module {
    // =========================================================================
    // DEBUG
    // Insert `tosa.identity` ops in the .tosa file to mark debug points.
    // This pattern rewrites them to tosa.custom with operator_name="debug_print"
    // so the runtime can print the tensor contents.
    //
    // Usage in .tosa:
    //   %dbg = tosa.identity %some_result : (tensor<?x28x1x40xi8>) -> tensor<?x28x1x40xi8>
    // =========================================================================

    pdl.pattern @debug_print : benefit(10) {
        %ty  = pdl.type
        %inp = pdl.operand : %ty

        %root = pdl.operation "tosa.identity"
                        (%inp : !pdl.value)
                        -> (%ty : !pdl.type)

        pdl.rewrite %root {
            %name  = pdl.attribute = "debug_print"
            %empty = pdl.attribute = ""
            %custom = pdl.operation "tosa.custom"
                                (%inp : !pdl.value)
                                {"operator_name"       = %name,
                                 "domain_name"          = %empty,
                                 "implementation_attrs" = %empty}
                                -> (%ty : !pdl.type)
            pdl.replace %root with %custom
        }
    }

//    // =========================================================================
//    // Layer 0: depthwise_conv2d
//    //   input:  tensor<?x30x1x40xi8>
//    //   weight: tensor<3x1x40x1xi8>   (HxWxICxM, H=3)
//    //   output: tensor<?x28x1x40xi32>
//    // =========================================================================
//
//    pdl.pattern @dw_conv2d_layer_0 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.depthwise_conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight shape 3x1x40x1: check H (dim 0) = 3
//        %h_dim  = pdl.attribute = 0 : i32
//        %h_size = pdl.attribute = 3 : i32
//        pdl.apply_native_constraint "check_type"(
//            %wTy, %h_dim, %h_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "dw_conv2d_layer_0"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 1: conv2d (pointwise 1x1)
//    //   input:  tensor<?x28x1x40xi8>
//    //   weight: tensor<128x1x1x40xi8>  (OCxHxWxIC, IC=40 — unique)
//    //   output: tensor<?x28x1x128xi32>
//    // =========================================================================
//
//    pdl.pattern @conv2d_layer_1 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight shape 128x1x1x40: check IC (dim 3) = 40
//        %ic_dim  = pdl.attribute = 3 : i32
//        %ic_size = pdl.attribute = 40 : i32
//        pdl.apply_native_constraint "check_type"(
//            %wTy, %ic_dim, %ic_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "conv2d_layer_1"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 2: depthwise_conv2d
//    //   input:  tensor<?x28x1x128xi8>
//    //   weight: tensor<5x1x128x1xi8>   (HxWxICxM, H=5)
//    //   output: tensor<?x24x1x128xi32>
//    // =========================================================================
//
//    pdl.pattern @dw_conv2d_layer_2 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.depthwise_conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight shape 5x1x128x1: check H (dim 0) = 5
//        %h_dim  = pdl.attribute = 0 : i32
//        %h_size = pdl.attribute = 5 : i32
//        pdl.apply_native_constraint "check_type"(
//            %wTy, %h_dim, %h_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "dw_conv2d_layer_2"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 3: conv2d (pointwise 1x1)
//    //   input:  tensor<?x24x1x128xi8>   <- H=24 distinguishes from layer 5
//    //   weight: tensor<128x1x1x128xi8>
//    //   output: tensor<?x24x1x128xi32>
//    // =========================================================================
//
//    pdl.pattern @conv2d_layer_3 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight is 128x1x1x128 — same as layer 5.
//        // Disambiguate by input spatial H (dim 1) = 24.
//        %h_dim  = pdl.attribute = 1 : i32
//        %h_size = pdl.attribute = 24 : i32
//        pdl.apply_native_constraint "check_type"(
//            %inpTy, %h_dim, %h_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "conv2d_layer_3"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 4: depthwise_conv2d
//    //   input:  tensor<?x24x1x128xi8>
//    //   weight: tensor<10x1x128x1xi8>   (HxWxICxM, H=10)
//    //   output: tensor<?x15x1x128xi32>
//    // =========================================================================
//
//    pdl.pattern @dw_conv2d_layer_4 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.depthwise_conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight shape 10x1x128x1: check H (dim 0) = 10
//        %h_dim  = pdl.attribute = 0 : i32
//        %h_size = pdl.attribute = 10 : i32
//        pdl.apply_native_constraint "check_type"(
//            %wTy, %h_dim, %h_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "dw_conv2d_layer_4"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 5: conv2d (pointwise 1x1)
//    //   input:  tensor<?x15x1x128xi8>   <- H=15 distinguishes from layer 3
//    //   weight: tensor<128x1x1x128xi8>
//    //   output: tensor<?x15x1x128xi32>
//    // =========================================================================
//
//    pdl.pattern @conv2d_layer_5 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight is 128x1x1x128 — same as layer 3.
//        // Disambiguate by input spatial H (dim 1) = 15.
//        %h_dim  = pdl.attribute = 1 : i32
//        %h_size = pdl.attribute = 15 : i32
//        pdl.apply_native_constraint "check_type"(
//            %inpTy, %h_dim, %h_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "conv2d_layer_5"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 6: depthwise_conv2d
//    //   input:  tensor<?x15x1x128xi8>
//    //   weight: tensor<15x1x128x1xi8>   (HxWxICxM, H=15)
//    //   output: tensor<?x1x1x128xi32>
//    // =========================================================================
//
//    pdl.pattern @dw_conv2d_layer_6 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.depthwise_conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight shape 15x1x128x1: check H (dim 0) = 15
//        %h_dim  = pdl.attribute = 0 : i32
//        %h_size = pdl.attribute = 15 : i32
//        pdl.apply_native_constraint "check_type"(
//            %wTy, %h_dim, %h_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "dw_conv2d_layer_6"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 7: conv2d (pointwise 1x1)
//    //   input:  tensor<?x1x1x128xi8>
//    //   weight: tensor<32x1x1x128xi8>   (OCxHxWxIC, OC=32 — unique)
//    //   output: tensor<?x1x1x32xi32>
//    // =========================================================================
//
//    pdl.pattern @conv2d_layer_7 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight shape 32x1x1x128: check OC (dim 0) = 32
//        %oc_dim  = pdl.attribute = 0 : i32
//        %oc_size = pdl.attribute = 32 : i32
//        pdl.apply_native_constraint "check_type"(
//            %wTy, %oc_dim, %oc_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "conv2d_layer_7"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
//
//    // =========================================================================
//    // Layer 8: conv2d (pointwise 1x1, final classification head)
//    //   input:  tensor<?x1x1x32xi8>
//    //   weight: tensor<3x1x1x32xi8>    (OCxHxWxIC, OC=3 — unique)
//    //   output: tensor<?x1x1x3xi32>
//    // =========================================================================
//
//    pdl.pattern @conv2d_layer_8 : benefit(1) {
//        %inpTy = pdl.type
//        %wTy   = pdl.type
//        %resTy = pdl.type
//
//        %inp       = pdl.operand : %inpTy
//        %w         = pdl.operand : %wTy
//        %bias      = pdl.operand
//        %input_zp  = pdl.operand
//        %weight_zp = pdl.operand
//
//        %stride   = pdl.attribute = array<i64: 1, 1>
//        %pad      = pdl.attribute = array<i64: 0, 0, 0, 0>
//        %dilation = pdl.attribute = array<i64: 1, 1>
//
//        %root = pdl.operation "tosa.conv2d"
//                        (%inp, %w, %bias, %input_zp, %weight_zp
//                            : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                        {"stride" = %stride, "pad" = %pad, "dilation" = %dilation}
//                        -> (%resTy : !pdl.type)
//
//        // Weight shape 3x1x1x32: check OC (dim 0) = 3
//        %oc_dim  = pdl.attribute = 0 : i32
//        %oc_size = pdl.attribute = 3 : i32
//        pdl.apply_native_constraint "check_type"(
//            %wTy, %oc_dim, %oc_size : !pdl.type, !pdl.attribute, !pdl.attribute
//        )
//
//        pdl.rewrite %root {
//            %customName = pdl.attribute = "conv2d_layer_8"
//            %empty      = pdl.attribute = ""
//            %custom = pdl.operation "tosa.custom"
//                                (%inp, %w, %bias, %input_zp, %weight_zp
//                                    : !pdl.value, !pdl.value, !pdl.value, !pdl.value, !pdl.value)
//                                {"operator_name"       = %customName,
//                                 "domain_name"          = %empty,
//                                 "implementation_attrs" = %empty}
//                                -> (%resTy : !pdl.type)
//            pdl.replace %root with %custom
//        }
//    }
}
