import os
import logging

logger = logging.getLogger('spam_application')
logger.setLevel(logging.DEBUG)


def create_folder(folder):
    logger.info('Creating Folder: ' + folder)
    if not os.path.exists(folder):
        os.mkdir(folder)
    else:
        logger.info('Folder is already existing!')
