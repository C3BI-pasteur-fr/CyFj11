test_that("ellipse gates follow plain geometry of FlowJo's vertices and rotation", {
  gate <- list(
    xVertices = list(47.2, 151), yVertices = list(67.9, 15.9), rotationAngle = 10.6,
    xAxis = list(parameterSpec = list(name = "FSC-A"),
                 transform = list(transformType = "Linear", minRange = 0, maxRange = 262144)),
    yAxis = list(parameterSpec = list(name = "SSC-A"),
                 transform = list(transformType = "Linear", minRange = 0, maxRange = 1525)))
  eg <- CyFj11:::convert_ellipse_gate(gate, "e", 0, -4000, schema_version = numeric_version("4.0.0"))

  set.seed(1)
  x <- runif(20000, 0, 262144); y <- runif(20000, 0, 1525)
  ff <- flowCore::flowFrame(cbind(`FSC-A` = x, `SSC-A` = y))
  got <- as.vector(flowCore::filter(ff, eg)@subSet)

  # plain geometry in 0-256 display space
  u <- x / 262144 * 256 - mean(c(47.2, 151)); v <- y / 1525 * 256 - mean(c(67.9, 15.9))
  a <- abs(151 - 47.2) / 2; b <- abs(15.9 - 67.9) / 2; th <- 10.6 * pi / 180
  want <- ((u * cos(th) + v * sin(th)) / a)^2 + ((-u * sin(th) + v * cos(th)) / b)^2 <= 1

  expect_gt(mean(got == want), 0.999)
  expect_gt(sum(want), 100)

  # older / unknown schema versions keep the legacy convention (backward compatible)
  legacy <- CyFj11:::convert_ellipse_gate(gate, "e", 0, -4000, schema_version = numeric_version("3.0.0"))
  unknown <- CyFj11:::convert_ellipse_gate(gate, "e", 0, -4000)
  expect_equal(legacy@cov, unknown@cov)
  expect_false(isTRUE(all.equal(legacy@cov, eg@cov)))
})

test_that("schema version reaches ellipse gates through extract_all_gates", {
  gate <- list(
    type = "ellipse",
    xVertices = list(47.2, 151), yVertices = list(67.9, 15.9), rotationAngle = 10.6,
    xAxis = list(parameterSpec = list(name = "FSC-A"),
                 transform = list(transformType = "Linear", minRange = 0, maxRange = 262144)),
    yAxis = list(parameterSpec = list(name = "SSC-A"),
                 transform = list(transformType = "Linear", minRange = 0, maxRange = 1525)))
  pd <- list(p1 = list(definition = list(name = list("Lym"), type = "gate",
                                         gateDefinition = gate, desyncTable = list())))
  new <- CyFj11:::extract_all_gates(pd, "s1", schema_version = numeric_version("4.0.0"))
  old <- CyFj11:::extract_all_gates(pd, "s1", schema_version = numeric_version("3.0.0"))
  expect_length(new, 1L)
  expect_length(old, 1L)
  expect_false(isTRUE(all.equal(new[[1]]@cov, old[[1]]@cov)))
})
