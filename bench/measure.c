// Run a command and report what it cost: `measure <report> <command> [arguments...]`
// writes the command's user seconds, system seconds, and peak resident set
// size in KB to the file <report>, and exits with the command's status.
//
// The biofast runner uses this instead of measuring from Python because
// Linux carries a process's peak resident set size across exec: a child of a
// Python process would never report less than Python's own footprint, which
// is larger than several of the programs measured. This parent is small.
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
	struct rusage usage;
	int status;
	pid_t child;
	FILE *report;

	if (argc < 3) {
		fprintf(stderr, "Usage: measure <report> <command> [arguments...]\n");
		return 2;
	}
	child = fork();
	if (child < 0) {
		perror("measure: fork");
		return 2;
	}
	if (child == 0) {
		execvp(argv[2], argv + 2);
		perror(argv[2]);
		_exit(127);
	}
	if (wait4(child, &status, 0, &usage) < 0) {
		perror("measure: wait4");
		return 2;
	}
	report = fopen(argv[1], "w");
	if (report == NULL) {
		perror(argv[1]);
		return 2;
	}
	fprintf(report, "%ld.%06ld %ld.%06ld %ld\n",
	        (long)usage.ru_utime.tv_sec, (long)usage.ru_utime.tv_usec,
	        (long)usage.ru_stime.tv_sec, (long)usage.ru_stime.tv_usec,
	        usage.ru_maxrss);
	fclose(report);
	return WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
}
