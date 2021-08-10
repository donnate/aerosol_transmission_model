library(tidyverse)
library(methods)
options(warn = -1) 

l2_norm = function(x, y) {
    return (sqrt(x^2 + y^2))
}

Person = setRefClass("Person", fields = c(infected="logical", dose="numeric"))

Person$methods(initialize = function(prevalence)
    {
        infected = runif(1) < prevalence
        initFields(infected=infected, dose=0)
    })

ExactModel = setRefClass("ExactModel", 
                         fields = c(n_servers="numeric", 
                                    distance_between_people="numeric", 
                                    distance_between_queues="numeric", 
                                    entry_rate="numeric", 
                                    service_rate="numeric", 
                                    prevalence="numeric", 
                                    elapsed_time="numeric", 
                                    duration="numeric",
                                    hazard_rates="list", 
                                    queues="list",
                                    exited_people="list"))

ExactModel$methods(initialize = function(n_servers, distance_between_people, 
                                         distance_between_queues, entry_rate, 
                                         service_time, duration, prevalence=0.01, 
                                         hazard_rate=0.00913) 
    {
        queues = vector(mode="list", length=n_servers)
        for (i in seq.int(n_servers)) {
            queues[[i]] = list()
        }

        exited_people = list()

        # Determine hazard rates at relevant distances from contagious individual
        hazard_rates = list()
        q_spread = floor(2 / distance_between_queues)
        p_spread = floor(2 / distance_between_people)
        for (di in seq.int(-q_spread, q_spread)) {
            for (dp in seq.int(-p_spread, p_spread)) {
                distance = l2_norm(di * distance_between_queues, dp * distance_between_people)
                if (distance < 2) {
                    hazard_rates = append(hazard_rates, list(c(di, dp, hazard_rate * exp(-log(2.02) * (distance - 0.5)))))
                }
            }
        }

        # Initialize the model
        initFields(n_servers=n_servers, 
                    distance_between_people=distance_between_people, 
                    distance_between_queues=distance_between_queues, 
                    entry_rate=entry_rate, 
                    service_rate=1/service_time, 
                    prevalence=prevalence, 
                    elapsed_time=0, 
                    duration=duration,
                    hazard_rates=hazard_rates, 
                    queues=queues,
                    exited_people=exited_people)
    })

ExactModel$methods(add_person = function() 
    {
        lengths = map(queues, length)
        shortest = which.min(lengths)
        queues[[shortest]] <<- append(queues[[shortest]], Person(prevalence))
    })

ExactModel$methods(remove_person = function(index)
    {
        exited_people <<- append(exited_people, queues[[index]][[1]])
        queues[[index]] <<- queues[[index]][-1]
    })

ExactModel$methods(infect = function(index, position, time_step)
    {
        for (vec in hazard_rates) {
            di = vec[[1]]
            dp = vec[[2]]
            hazard_rate = vec[[3]]
            if (index + di > 0 & index + di <= length(queues)) {
                if (position + dp > 0 & position + dp <= length(queues[[index + di]])) {
                    queues[[index + di]][[position + dp]]$dose <<- queues[[index + di]][[position + dp]]$dose + hazard_rate * time_step
                }
            }
            # try(queues[[index + di]][[position + dp]]$dose <<- queues[[index + di]][[position + dp]]$dose + hazard_rate * time_step, silent=TRUE)
        }
        # for (idx in seq_along(queues)) {
        #     for (pos in seq_along(queues[[idx]])) {
        #         distance = l2_norm((index - idx) * distance_between_queues,
        #                             (position - pos) * distance_between_people)
        #         if (distance < 2) {
        #             queues[[idx]][[pos]]$dose <<- queues[[idx]][[pos]]$dose + 0.00913 * time_step * exp(-log(2.02) * (distance - 0.5))
        #         }
        #     }
        # }
    })

ExactModel$methods(step = function()
    {
        # Generate time to next event.
        busy_queues = which(sapply(queues, function(l) length(l) != 0))
        total_rate = length(busy_queues) * service_rate + entry_rate
        time_step = rexp(1, total_rate)

        # Check to see if event duration has been reached.
        if (elapsed_time + time_step > duration) {
            time_step = duration - elapsed_time
            elapsed_time <<- duration
        
            # Infect susceptible individuals.
            for (index in seq_along(queues)) {
                for (position in seq_along(queues[[index]])) {
                    if (queues[[index]][[position]]$infected) {
                        infect(index, position, time_step)
                    }
                }
            }

            return (FALSE)
        }

        else {
            elapsed_time <<- elapsed_time + time_step

            # Infect susceptible individuals.
            for (index in seq_along(queues)) {
                for (position in seq_along(queues[[index]])) {
                    if (queues[[index]][[position]]$infected) {
                        infect(index, position, time_step)
                    }
                }
            }

            # Determine whether the event is an entry...
            p_entry = entry_rate / total_rate
            if (runif(1) < p_entry) {
                add_person()
            }
            # ... or an exit.
            else {
                index_serviced = busy_queues[sample(1:length(busy_queues), 1)]
                remove_person(index_serviced)
            }

            return (TRUE)
        }
    })

ExactModel$methods(new_infections = function()
    {
        cumulative_sum = 0
        for (queue in queues) {
            for (person in queue) {
                if (!person$infected) {
                    cumulative_sum = cumulative_sum + 1 - exp(-person$dose)
                }
            }
        }

        for (person in exited_people) {
            if (!person$infected) {
                cumulative_sum = cumulative_sum + 1 - exp(-person$dose)
            }
        }

        return (cumulative_sum)
    })

simulation = function(n_trials, ...) {
    data = c()
    for (i in seq.int(n_trials)) {
        print(i)
        m = ExactModel(...)
        while (m$step()) {}
        n_data = m$new_infections()
        data = append(data, m$new_infections())
    }
    return (data)
}