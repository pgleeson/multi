job_id=$1


while IFS= read -r line; do

    if [[ $line == *"="* ]]
    then
      lineclean=$( echo ${line} | xargs)
      var=$(echo $lineclean | cut -d'=' -f 1 )
      val=$(echo $lineclean | cut -d'=' -f 2 )
      #echo "Setting: $var to $val ...";
      export $var=$val
    fi
    
done < ~/nsgrest.admin.osb

echo "  Read in data for NSG REST API user: $ADMIN_USERNAME to access: $URL..."

url=$URL/job/$UMBRELLA_APPNAME.$USER_USERNAME/$job_id

echo "Accessing: $url for user $ADMIN_USERNAME on behalf of $USER_USERNAME"

xml_str=$(curl -u $ADMIN_USERNAME:$ADMIN_PASSWORD \
     -H cipres-appkey:$UMBRELLA_APPID \
     -H cipres-eu:$USER_USERNAME \
     -H cipres-eu-email:$USER_EMAIL \
     -H cipres-eu-institution:$USER_INSTITUTION \
     -H cipres-eu-country:$USER_COUNTRY \
      $url)
      
__='
echo "==============="
echo "[$xml_str]"
echo "==============="
'

job_stage=$(xmllint --xpath "string(//jobStage)"  - <<<"$xml_str")
echo "Done checking job: " $job_stage

if [ $job_stage = "COMPLETED" ]
then
    echo "Downloading results of " $job_id

    xml_str=$(curl -u $ADMIN_USERNAME:$ADMIN_PASSWORD \
         -H cipres-appkey:$UMBRELLA_APPID \
         -H cipres-eu:$USER_USERNAME \
         -H cipres-eu-email:$USER_EMAIL \
         -H cipres-eu-institution:$USER_INSTITUTION \
         -H cipres-eu-country:$USER_COUNTRY \
          $URL/job/$UMBRELLA_APPNAME.$USER_USERNAME/$job_id/output)
          

    output_url=$(xmllint --xpath "string(/results/jobfiles/jobfile[parameterName='outputfile']/downloadUri/url)"  - <<<"$xml_str")

    echo "Results are in: "$output_url

    mv output.tar.gz /tmp

    curl -u $ADMIN_USERNAME:$ADMIN_PASSWORD \
         -H cipres-appkey:$UMBRELLA_APPID \
         -H cipres-eu:$USER_USERNAME \
         -H cipres-eu-email:$USER_EMAIL \
         -H cipres-eu-institution:$USER_INSTITUTION \
         -H cipres-eu-country:$USER_COUNTRY \
          -O -J $output_url
          
    mkdir $job_id
    mv output.tar.gz $job_id
    cd $job_id
    tar xvzf output.tar.gz
    
    echo
    echo "Results should now be in ./$job_id/temp"
    echo
    
else

    echo
    echo "Job" $job_id "is not finished yet!"
    echo
fi
