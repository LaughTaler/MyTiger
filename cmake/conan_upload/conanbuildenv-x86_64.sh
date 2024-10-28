script_folder="/home/yhf/code/tiger1.5/testUploadConan"
echo "echo Restoring environment" > "$script_folder/deactivate_conanbuildenv-x86_64.sh"
for v in 
do
    is_defined="true"
    value=$(printenv $v) || is_defined="" || true
    if [ -n "$value" ] || [ -n "$is_defined" ]
    then
        echo export "$v='$value'" >> "$script_folder/deactivate_conanbuildenv-x86_64.sh"
    else
        echo unset $v >> "$script_folder/deactivate_conanbuildenv-x86_64.sh"
    fi
done

