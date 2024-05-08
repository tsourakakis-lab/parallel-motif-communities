licenses(["notice"])

cc_library(
    name = "Thresholding",
    hdrs = ["Thresholding.h"],
    deps = [
    	"//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
        "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
        "//gbbs",
    ],
)


cc_binary(
    name = "Thresholding_main",
    srcs = ["Thresholding.cc"],
    deps = [
        ":Thresholding"
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
