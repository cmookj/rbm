#
# Dependencies
#
include(FetchContent)

# Google Test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG 52eb8108c5bdec04579160ae17225d66034bd723 # release-1.17.0
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# vma
FetchContent_Declare(
  vma
  GIT_REPOSITORY https://github.com/cmookj/vma.git
  GIT_TAG origin/main
)
FetchContent_MakeAvailable(vma)

# toolbox
FetchContent_Declare(
  toolbox
  GIT_REPOSITORY https://github.com/cmookj/toolbox.git
  GIT_TAG origin/main
)
FetchContent_MakeAvailable(toolbox)
