# Update ecoli_fbi submodule
echo "Updating ecoli_fbi submodule to the most current commit"
git submodule update --init --recursive
pushd bifrost_sp_ecoli/ecoli_fbi || exit

# Fetch the latest changes from the remote repository
git fetch origin --tags
git checkout main  # or the specific branch you want to track
git pull origin main --tags

# Show the new commit hash after the update
LATEST_COMMIT=$(git rev-parse HEAD)
echo "Updated commit hash of ecoli_fbi after update: $LATEST_COMMIT"

# print the latest tag of ecoli_fbi
LATEST_TAG_COMMIT=$(git for-each-ref --sort=-creatordate refs/tags|head -1|cut -f1 -d ' ')
LATEST_TAG=$(git for-each-ref --sort=-creatordate refs/tags|head -1|cut -f3 -d '/')

echo "Checking commit hash for the latest tag of ecoli_fbi: $LATEST_TAG_COMMIT"
echo "Checking the lastest tag of ecoli_fbi: $LATEST_TAG"
  
git checkout $LATEST_TAG_COMMIT

# Navigate back to the main project directory
popd
