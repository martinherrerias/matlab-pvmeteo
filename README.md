# matlab-pvmeteo
Functions for manipulation of meteorological data in solar applications.

## Use
#### 1. Add all folders to the matlab path:

`addpath(genpath('./matlab-pvmeteo'));`

#### 2. Add user accounts:

User accounts are required to access [SoDa](http://www.soda-pro.com/) and [geonames](https://www.geonames.org/). In most cases these provide optional features (look for the option `useweb = false`). To use the services, modify the file `./Remote/useraccounts.txt` with your credentials and save as `useraccounts.m`:

```
copyfile('./Remote/useraccounts.txt','./Remote/useraccounts.m');
% (update `useraccounts.m` with your credentials)
```

To avoid repeated web requests (or to allow manual downloads), _mcclear*.csv_ and _merra2*.csv_ files can also be retrieved from local folders in the search-path. Files are identified by name (see `getremote.searchexisting`). The functions `merra2_healthcheck` and `mcclear_healthcheck` can be used periodically to auto-rename the files and speed-up the search:

```
addpath('~/Documents/CAMS/MERRA2');
copyfile('./Remote/resources/merra2_healthcheck.m','~/Documents/CAMS/MERRA2');

addpath('~/Documents/CAMS/McClear');
copyfile('./Remote/resources/mcclear_healthcheck.m','~/Documents/CAMS/McClear');
```
