from bioflow.utils.general_utils.high_level_os_io import mkdir_recursive
from bioflow.utils.general_utils.internet_io import url_to_local, marbach_post_proc
from bioflow.utils.log_behavior import get_logger
from bioflow.configs.main_configs import DB_locations, sources_location


log = get_logger(__name__)


def pull_online_dbs(_db_locations: dict =DB_locations,
                    _sources_location: dict =sources_location) -> None:
    """
    Pulls the databases mentionned online to the directories specified in the configs file

    :param _db_locations:
    :param sources_location:
    :return:
    """

    for DB_type, location_dict in _db_locations.items():

        # print(DB_type, location_dict)

        if 'inactive' in list(location_dict.keys()) and location_dict['inactive'] == 'True':
            continue

        local = location_dict['local'].replace('$DB_HOME$', _sources_location)
        onlines = [location.strip() for location in location_dict['online'].split(',')]

        if 'rename' in list(location_dict.keys()):
            renames = [location.strip() for location in location_dict['rename'].split(',')]

            for online, rename in zip(onlines, renames):
                log.info('loading %s database from %s to %s',
                         DB_type.lower(), location_dict['online'], local)
                mkdir_recursive(local)
                url_to_local(online, local, rename=rename)

        else:
            for online in onlines:
                log.info('loading %s database from %s to %s',
                         DB_type.lower(), location_dict['online'], local)
                mkdir_recursive(local)
                url_to_local(online, local)

        if DB_type == 'MARBACH':
            log.info('cleaning up local marbach databse')
            marbach_post_proc(local)
