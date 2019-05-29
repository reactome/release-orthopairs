
releaseNumber=
configFilepath=src/main/resources/config.properties

while (( "$#" )); do
	case "$1" in
		-r|--release)
			releaseNumber=$2;
			shift 2;
			;;
		-c|--config)
			configFilepath=$2;
			shift 2;
			;;
		-*|--*=)
			echo "Error: Unsupported flag $1";
			exit 1;
	esac
done
echo $releaseNumber;
if [ -z "$releaseNumber" ]
then
	echo "Updates config.properties with release information used by Orthopairs"
	echo "Usage: bash updateOrthopairsConfig.sh -r releaseNumber [-c configFilepath]"
	exit 1
fi

sed -i "s/releaseNumber=/releaseNumber=${releaseNumber}/g" $configFilepath
