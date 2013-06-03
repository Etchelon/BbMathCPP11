#include "BbStopWatch.hpp"
#include <utility>

using namespace std;

void BbStopWatch::start()
{
	if (!isRunning)
	{
		splitTime = startTime = hr_clock::now();
		isRunning = true;
	}
}

void BbStopWatch::stop()
{
	if (isRunning)
	{
		stopTime = hr_clock::now();
		isRunning = false;
	}
}

void BbStopWatch::reset()
{
	stopTime = splitTime = startTime = hr_clock::now();
	isRunning = false;
}

void BbStopWatch::restart()
{
	reset();
	start();
}

double BbStopWatch::get_duration(const hr_t_point& e, const hr_t_point& s, const TimeUOM& p) const
{
	using namespace chrono;

	switch (p)
	{
	case TimeUOM::hours:
		return duration_cast<duration<double, std::ratio<3600>>>(e - s).count();

	case TimeUOM::minutes:
		return duration_cast<duration<double, std::ratio<60>>>(e - s).count();

	case TimeUOM::seconds:
		return duration_cast<duration<double, std::ratio<1>>>(e - s).count();

	case TimeUOM::milliseconds:
		return duration_cast<duration<double, std::milli>>(e - s).count();

	case TimeUOM::microseconds:
		return duration_cast<duration<double, std::micro>>(e - s).count();

	case TimeUOM::nanoseconds:
		return duration_cast<duration<double, std::nano>>(e - s).count();
	}
}

double BbStopWatch::split()
{
	if (!isRunning)
	{
		return 0.;
	}

	auto now = hr_clock::now();
	double duration = get_duration(now, splitTime, period);

	splitTime = std::move(now);

	return duration;
}

double BbStopWatch::get_total_running_time() const
{
	auto now = isRunning ? hr_clock::now() : stopTime;
	double duration = get_duration(now, startTime, period);

	return duration;
}
