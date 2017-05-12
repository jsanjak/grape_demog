/*!
  \brief Testing ground for more flexible API to evolve populations -- single
  deme version.
*/
#ifndef __CLONAL_SAMPLE_DIPLOID_CLONAL_HPP__
#define __CLONAL_SAMPLE_DIPLOID_CLONAL_HPP__

#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/dispatch.hpp>
#include <fwdpy11/rules/wf_rules.hpp>
namespace clonal
{
        // single deme, N changing
        //! \brief Experimental variant where the population rules are
        //! implemented via an external policy
        template <typename gamete_type, typename gcont_t_allocator,
                  typename mutation_type, typename mcont_t_allocator,
                  typename diploid_geno_t,
                  typename diploid_vector_type_allocator,
                  typename diploid_fitness_function, typename mutation_model,
                  typename recombination_policy,
                  template <typename, typename> class gcont_t,
                  template <typename, typename> class mcont_t,
                  template <typename, typename> class diploid_vector_type,
                  typename popmodel_rules = fwdpy11::wf_rules,
                  typename mutation_removal_policy = std::true_type,
                  typename gamete_insertion_policy = KTfwd::emplace_back>
        double
        sample_diploid_clonal(
            const gsl_rng *r, gcont_t<gamete_type, gcont_t_allocator> &gametes,
            diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
                &diploids,
            mcont_t<mutation_type, mcont_t_allocator> &mutations,
            std::vector<KTfwd::uint_t> &mcounts, const unsigned &N_curr,
            const unsigned &N_next, const double &mu,
            const mutation_model &mmodel, const recombination_policy &rec_pol,
            const diploid_fitness_function &ff,
            typename gamete_type::mutation_container &neutral,
            typename gamete_type::mutation_container &selected,
            const double &f = 0., popmodel_rules &&pmr = popmodel_rules(),
            const mutation_removal_policy &mp = mutation_removal_policy(),
            const gamete_insertion_policy &gpolicy_mut
            = gamete_insertion_policy())
        {
            assert(N_curr == diploids.size());

            auto gamete_recycling_bin
                = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
            auto mutation_recycling_bin
                = KTfwd::fwdpp_internal::make_mut_queue(mcounts);

            KTfwd::experimental::dispatch_w(pmr, diploids, gametes, mutations, ff);

#ifndef NDEBUG
            for (const auto &g : gametes)
                assert(!g.n);
#endif
            const auto parents(diploids); // copy the parents

            // Change the population size
            if (diploids.size() != N_next)
                {
                    diploids.resize(N_next);
                }
            for (auto &dip : diploids)
                {
                    // Pick parent 1
                    size_t p1 = pmr.pick1(r);
                    assert(p1 < parents.size());

                    std::size_t p1g1 = parents[p1].first;
                    std::size_t p1g2 = parents[p1].second;

                    dip.first = p1g1; 
                    dip.second = p1g2; 

                    gametes[dip.first].n++;
                    gametes[dip.second].n++;

                    // now, add new mutations
                    dip.first = KTfwd::mutate_gamete_recycle(
                        mutation_recycling_bin, gamete_recycling_bin, r, mu,
                        gametes, mutations, dip.first, mmodel, gpolicy_mut);
                    dip.second = KTfwd::mutate_gamete_recycle(
                        mutation_recycling_bin, gamete_recycling_bin, r, mu,
                        gametes, mutations, dip.second, mmodel, gpolicy_mut);

                    assert(gametes[dip.first].n);
                    assert(gametes[dip.second].n);
                    KTfwd::experimental::dispatch_update(pmr, r, dip, parents[p1], parents[p1],
                                    gametes, mutations, ff);
                }
#ifndef NDEBUG
            for (const auto &dip : diploids)
                {
                    assert(gametes[dip.first].n);
                    assert(gametes[dip.first].n <= 2 * N_next);
                    assert(gametes[dip.second].n);
                    assert(gametes[dip.second].n <= 2 * N_next);
                }
#endif
            KTfwd::fwdpp_internal::process_gametes(gametes, mutations, mcounts);
            assert(KTfwd::popdata_sane(diploids, gametes, mutations, mcounts));
            KTfwd::fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts,
                                           2 * N_next, mp);
            assert(KTfwd::check_sum(gametes, 2 * N_next));
            return pmr.wbar;
        }

        // single deme, N constant

        //! \brief Experimental variant where the population rules are
        //! implemented via an external policy
        template <typename gamete_type, typename gcont_t_allocator,
                  typename mutation_type, typename mcont_t_allocator,
                  typename diploid_geno_t,
                  typename diploid_vector_type_allocator,
                  typename diploid_fitness_function, typename mutation_model,
                  typename recombination_policy,
                  template <typename, typename> class gcont_t,
                  template <typename, typename> class mcont_t,
                  template <typename, typename> class diploid_vector_type,
                  typename popmodel_rules = fwdpy11::wf_rules ,
                  typename mutation_removal_policy = std::true_type,
                  typename gamete_insertion_policy = KTfwd::emplace_back>
        double
        sample_diploid_clonal(
            const gsl_rng *r, gcont_t<gamete_type, gcont_t_allocator> &gametes,
            diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
                &diploids,
            mcont_t<mutation_type, mcont_t_allocator> &mutations,
            std::vector<KTfwd::uint_t> &mcounts, const unsigned &N_curr,
            const double &mu, const mutation_model &mmodel,
            const recombination_policy &rec_pol,
            const diploid_fitness_function &ff,
            typename gamete_type::mutation_container &neutral,
            typename gamete_type::mutation_container &selected,
            const double &f = 0., popmodel_rules &&pmr = popmodel_rules(),
            const mutation_removal_policy &mp = mutation_removal_policy(),
            const gamete_insertion_policy &gpolicy_mut
            = gamete_insertion_policy())
        {
            return sample_diploid_clonal(
                r, gametes, diploids, mutations, mcounts, N_curr, N_curr, mu,
                mmodel, rec_pol, ff, neutral, selected, f, pmr, mp,
                gpolicy_mut);
        }
}
#endif
