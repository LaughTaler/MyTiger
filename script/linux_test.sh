#!/bin/bash

# 调用程序
# cp /home/gitlab-runner/release/AFT_Tri_test /home/gitlab-runner/tiger/aft_tri_test/
cd /home/gitlab-runner/tiger/aft_tri_test
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/gitlab-runner/tiger/release
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/gitlab-runner/tiger/aft_tri_test
# ./AFT_Tri_test [test1]

# # 获取返回值
# return_value=$?

# # 根据返回值进行处理
# if [ $return_value -eq 0 ]; then
#     echo "Program executed successfully."
# else
#     echo "Program failed with return value: $return_value"
#     exit $return_value
# fi

# 调用程序
./AFT_Tri_test [test2]

# 获取返回值
return_value=$?

# 根据返回值进行处理
if [ $return_value -eq 0 ]; then
    echo "Program executed successfully."
else
    echo "Program failed with return value: $return_value"
    exit $return_value
fi
# 调用程序
./AFT_Tri_test [test3]

# 获取返回值
return_value=$?

# 根据返回值进行处理
if [ $return_value -eq 0 ]; then
    echo "Program executed successfully."
else
    echo "Program failed with return value: $return_value"
    exit $return_value
fi
# 调用程序
./AFT_Tri_test [test4]

# 获取返回值
return_value=$?

# 根据返回值进行处理
if [ $return_value -eq 0 ]; then
    echo "Program executed successfully."
else
    echo "Program failed with return value: $return_value"
    exit $return_value
fi

# cd ../dt_remesh_test

# ./DT_Remesh_test [test]

# # 获取返回值
# return_value=$?

# # 根据返回值进行处理
# if [ $return_value -eq 0 ]; then
#     echo "Program executed successfully."
# else
#     echo "Program failed with return value: $return_value"
#     exit $return_value
# fi

exit $return_value