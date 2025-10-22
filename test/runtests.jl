using TestItemRunner

@run_package_tests verbose = true filter = test -> (:skipci ∉ test.tags)
