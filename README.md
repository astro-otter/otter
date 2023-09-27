# OTTER
### **O**pen **T**idal disrup**T**ion **E**vent **R**epository

## Repo Organization
| Directory | Contents |
|------------|------------|
| `data` | Testing data for the catalog |
| `db` | scripts to create the arangdb database |
| `py/otter` | A pip installable API for interfacing with the ArangoDB database|
| `webotter` | The flask app directory |

## Developer Instructions
1. Clone this repo:
   ```
   git clone https://github.com/noahfranz13/otter.git
   ```	 
2. Build the database. First install arangodb from
   https://www.arangodb.com/download-major/ubuntu/.
   Then, you can build the database by running the
   following commands:
   ```
   cd db
   ./setup.sh
   ```
3. Install tidecat, the API for this database. From
   the root directory of this repo run
   ```
   pip install -e .
   ```
4. Open the website locally. First you have to navigate
   into the directory hosting the flask info and then
   you can build the website. Run the following commands
   ```
   cd webtide
   flask --app webtide run
   ```

