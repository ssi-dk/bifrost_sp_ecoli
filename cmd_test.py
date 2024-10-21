
import pymongo
from pymongo.mongo_client import MongoClient
from pymongo.server_api import ServerApi
import os
import hashlib

def pingclient(key):
    """
    Connect to the MongoDB server and ping it to confirm the connection.

    Parameters:
    key (str): The MongoDB URI to connect to the server.

    Raises:
    ValueError: Key is None, indicating that the MongoDB URI is not set.
    """
    
    if key is None:
        raise ValueError("MongoDB URI is not set in the environment variable 'BIFROST_DB_KEY'. ")


    # Create a new client and connect to the server
    client = MongoClient(key,ssl=True)
    print("client work")
    # Send a ping to confirm a successful connection
    try:
        client.admin.command('ping')
        print("Pinged your deployment. You successfully connected to MongoDB!")
    except Exception as e:
        print(e)

def calculate_md5(data):
    """
    Calculate the MD5 checksum of the provided data.

    Parameters:
    data (str | list): The input data to calculate the MD5 checksum.

    Returns:
    str: The MD5 checksum of the input data as a hexadecimal string.

    Raises:
    ValueError: If the input data is neither a string nor a list.
    """

    if isinstance(data, list):
        # Convert the list to a string
        data_str = str(data)
    elif isinstance(data, str):
        data_str = data
    else:
        raise ValueError("Input must be a string or a list.")
    
    # Calculate and return the MD5 checksum
    return hashlib.md5(data_str.encode()).hexdigest()

def database_names(mongo_client,printcol=0,outputfile=None,outputtype=None):
    """
    Retrieve and optionally print the names of all databases from the MongoDB client.

    Parameters:
    mongo_client (MongoClient): The MongoDB client to retrieve names from.
    printcol (int): Flag to control printing of database details (1 to print, 0 to skip).
    outputfile (str): Optional; file name to write the results.
    outputtype (str): Optional; file mode ('w' for write, 'a' for append).

    Returns:
    list: A list containing the names of the databases.
    """

    database_names = mongo_client.list_database_names()

    i = 0
    for i, db_name in enumerate(database_names):
        if printcol == 1:
            print(f"Database no. {i} with name: {db_name}")

        # Store MD5 of collections list
    md5_db_names = calculate_md5(str(database_names))
    
    if printcol == 1:
        print(f"MD5 before: {md5_db_names}")
        print(f"Collections names: {database_names}")

    if outputfile:
        if outputtype:    
            with open(outputfile, outputtype) as f:
                f.write(f"{md5_db_names}\tDatabases:\t{database_names}\n")  # Write MD5 checksum as the first colum 
        else:
            print("specify file mode, w: write, a: append")

    return database_names

def collection_names(mongo_client,database,printcol=0,outputfile=None,outputtype=None):
    """
    Retrieve and optionally print the names of all collections in a specified database.

    Parameters:
    mongo_client (MongoClient): The MongoDB client to retrieve collections from.
    database (str): The name of the database to list collections from.
    printcol (int): Flag to control printing of database details (1 to print, 0 to skip).
    outputfile (str): Optional; file name to write the results.
    outputtype (str): Optional; file mode ('w' for write, 'a' for append).

    Returns:
    list: A list containing the names of the collections in the specified database.
    """
    db = mongo_client.get_database(database)
    collections_names = db.list_collection_names()
    
    # Store MD5 of collections list
    md5_names = calculate_md5(str(collections_names))
    
    if printcol == 1:
        print(f"MD5 before: {md5_names}")
        print(f"Collections names: {collections_names}")

    if outputfile:
        if outputtype:    
            with open(outputfile, outputtype) as f:
                f.write(f"{md5_names}\tCollection:\t{collections_names}\n")  # Write MD5 checksum as the first colum 
        else:
            print("specify file mode, w: write, a: append")

    return collections_names

def entry_names(mongo_client,database,printcol=0,outputfile=None,outputtype=None):
    """
    Retrieve and optionally print the names of all entries in each collection of a specified database.

    Parameters:
    mongo_client (MongoClient): The MongoDB client to retrieve entries from.
    database (str): The name of the database to list entries from.
    printcol (int): Flag to control printing of entry details (1 to print, 0 to skip).
    outputfile (str): Optional; file name to write the results.
    outputtype (str): Optional; file mode ('w' for write, 'a' for append).

    Returns:
    list: A list containing the names of the entries (keys) from the first documents in each collection.
    """

    db = mongo_client.get_database(database)
    collections_before = db.list_collection_names()

    entriesnames=[]
    for collection_name in collections_before:
        collection = db[collection_name]
        first_entry = collection.find_one()  # Retrieves the first document
        if printcol == 1:
            print("collection ",collection)
            if first_entry:
                print(f"First entry in '{collection_name}' collection: {first_entry}")
            else:
                print(f"No entries found in '{collection_name}' collection.")
        
        for entries in first_entry.keys():
            if printcol == 1:
                print("entries ",entries)
            entriesnames.append(entries)

    md5_entries = calculate_md5(entriesnames)

    colnames = list(first_entry.keys())
    if outputfile:
        if outputtype:    
            with open(outputfile, outputtype) as f:
                f.write(f"{md5_entries}\tEntries:\t{colnames}\t\n")  # Write MD5 checksum as the first colum 
        else:
            print("specify file mode, w: write, a: append")

    return entriesnames

if __name__ == "__main__":
    uri = os.environ["BIFROST_DB_KEY"] 

    print("Key works \n")
    print(uri)
    pingclient(uri)

    client = pymongo.MongoClient(uri)