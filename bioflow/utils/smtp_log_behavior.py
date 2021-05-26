"""
This modules carries the smtp logging logic
"""
import logging
from logging.handlers import SMTPHandler
from smtplib import SMTP
from datetime import datetime
from email.message import EmailMessage
import platform
import atexit

from bioflow.configs.main_configs import smtp_logging_parameters, smtp_logging

if not smtp_logging:
    raise Exception('Please configure the smtp logging module before trying to use it')

# define a formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s')

local_host = smtp_logging_parameters['local_host']
local_mail_account = smtp_logging_parameters['local_mail_account']
reporting_target_mail = smtp_logging_parameters['reporting_target_mail']

# debug
print('smtp logging enabled. '
      'Will be using local host %s with local mail address %s to report to %s.' %
      (local_host, local_mail_account, reporting_target_mail,))

input('press entre to continue')

mail_handler = SMTPHandler(mailhost=local_host,
                           fromaddr=local_mail_account,
                           toaddrs=reporting_target_mail,
                           subject="BioFlow runtime error")

mail_handler.setLevel(logging.CRITICAL)

start_date = datetime.now()
current_device = platform.node()

global _success
_success = True

def started_process():
    """
    Reports the start of a run

    :return: None
    """
    # print(smtp_logging_parameters['local_host'])

    mime_message = EmailMessage()
    mime_message['From'] = local_mail_account
    mime_message['To'] = reporting_target_mail

    with SMTP(host=local_host) as smtp_server:
        mime_message['Subject'] = 'Run started: BioFlow on %s on %s' % (current_device,
                                                               start_date.strftime('%x %X'))
        mime_message.set_content('Run has started successfully')

        smtp_server.send_message(mime_message)


def completed():
    """
    Reports the successful completion of a run

    :return: None
    """

    mime_message = EmailMessage()
    mime_message['From'] = local_mail_account
    mime_message['To'] = reporting_target_mail

    with SMTP(host=local_host) as smtp_server:
        mime_message['Subject'] = 'Run completed: BioFlow on %s on %s' % (current_device,
                                                               start_date.strftime('%x %X'))
        mime_message.set_content('Run has completed after %s minutes' %
                                 ((datetime.now() - start_date).total_seconds()//60.))

        smtp_server.send_message(mime_message)


def smtp_error_bail_out():
    """
    Reports an error of the SMTPHandler when the handler itself fails.

    :return:
    """

    mime_message = EmailMessage()
    mime_message['From'] = local_mail_account
    mime_message['To'] = reporting_target_mail

    with SMTP(host=local_host) as smtp_server:
        mime_message['Subject'] = 'SMTPHandler error bail-out on %s' % current_device
        mime_message.set_content("There was an error in the SMTP handler code at for run "
                                 "started on %s on machine %s resulting in a bail-out."
                                 % (datetime.now().strftime('%x %X'), current_device))
        # print(smtp_server.noop())
        smtp_server.send_message(mime_message)



# def smtp_logger(message):  # wrap the logger into a try/except loop
#     try:
#         _raw_smtp_logger.critical(message)
#     except Exception as e:
#         smtp_error_bail_out()
#         raise e

atexit.register(completed)  # and register the run termination report

started_process()

# will interfere if we now have the logger. When the smtp is enabled we will always be
# raising logging exception handler
