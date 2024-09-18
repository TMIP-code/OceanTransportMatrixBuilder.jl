using TestItemRunner

@run_package_tests verbose=true filter=ti->!(:skipci in ti.tags)
