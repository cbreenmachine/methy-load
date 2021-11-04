
echo "Her's some message. Continue? y/n"
read -n 1 k <&1
if [[ $k = n ]] ;
then
exit ERRCODE "Re-configure directories or delete files as needed!"
elif [[ $k = y ]]; then
echo "continuing"

fi