set -e

# Build the NML

python ISN.py -test

# Generate Neuron or NetPyNE in temp folder
cd temp
jnml LEMS_ISN_net.xml -netpyne

# create an init.py
cp LEMS_ISN_net_netpyne.py init.py




######################################
#   The rest here should not need to be altered, as long as everything (incl init.py) is in temp/

cd ..

zip -r input temp/*

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


url=$URL/job/$UMBRELLA_APPNAME.$USER_USERNAME

echo "Accessing: $url for user $ADMIN_USERNAME on behalf of $USER_USERNAME"

xml_str=$(curl -u $ADMIN_USERNAME:$ADMIN_PASSWORD \
     -H cipres-appkey:$UMBRELLA_APPID \
     -H cipres-eu:$USER_USERNAME \
     -H cipres-eu-email:$USER_EMAIL \
     -H cipres-eu-institution:$USER_INSTITUTION \
     -H cipres-eu-country:$USER_COUNTRY \
     -F tool='OSBPYNEURON74' \
     -F input.infile_=@input.zip \
     -F metadata.clientJobId=ISN \
     -F metadata.statusEmail=true \
     -F vparam.number_cores_=8 \
     -F vparam.number_nodes_=1 \
      $url)
      
      
job_id=$(xmllint --xpath "string(//jobHandle)" - <<<"$xml_str")
echo "Done submitting job: " $job_id
echo ""
echo "Retrieve results with:"
echo
echo "    ./retrieve_from_nsg.sh" $job_id
echo
      


