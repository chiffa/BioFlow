__author__ = 'ank'

import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='Go_UP_insert_log.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

def clean_with_log(ObjectType):
    i=0
    for elt in ObjectType.get_all():
        ID=str(elt).split('/')[-1][:-1]
        i+=1
        logging.debug('del %s', i)
        ObjectType.delete(ID)

def clean(ObjectType):
    i=0
    for elt in ObjectType.get_all():
        ID=str(elt).split('/')[-1][:-1]
        i+=1
        ObjectType.delete(ID)
