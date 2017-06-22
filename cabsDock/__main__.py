from job import Job


def run_job():
    job_args = {}
    j = Job(**job_args)
    j.run_job()

if __name__ == '__main__':
    run_job()
