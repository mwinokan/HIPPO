#!/bin/bash

### RESETS CACHES, MIGRATIONS, AND HIPPO DATABASE

rm -r ../web/__pycache__
rm -r ../hippo/__pycache__
rm -r ../hippo/migrations
rm HIPPO-DB.sqlite3
