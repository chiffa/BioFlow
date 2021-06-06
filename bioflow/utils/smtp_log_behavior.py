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

local_host = smtp_logging_parameters['local_host']
local_mail_account = smtp_logging_parameters['local_mail_account']
reporting_target_mail = smtp_logging_parameters['reporting_target_mail']

# define a formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s')

start_date = datetime.now()
current_device = platform.node()

report_name = 'Run critical error: BioFlow on %s on %s' % (current_device,
                                                          start_date.strftime('%x %X'))

mail_handler = SMTPHandler(mailhost=local_host,
                           fromaddr=local_mail_account,
                           toaddrs=reporting_target_mail,
                           subject=report_name)

mail_handler.setLevel(logging.CRITICAL)


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
        smtp_server.send_message(mime_message)


atexit.register(completed)  # and register the run termination report

started_process()

