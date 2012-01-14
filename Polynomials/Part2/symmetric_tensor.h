#ifndef SYMMETRIC_TENSOR_H
#define SYMMETRIC_TENSOR_H

#include <array>
#include <utility>
#include <iostream>

//Work around for broken const expr on g++ 4.4
template<int n>
struct Factorial { enum { value = n * Factorial<n-1>::value }; };

template<>
struct Factorial<0> { enum { value = 1 }; };

template<int n, int r>
struct Binomial { enum { value = Factorial<n>::value / (Factorial<r>::value * Factorial<n-r>::value) }; };

//Combinatorial helper methods
inline int factorial(int n) {
    int r = 1;
    for(int i=2; i<=n; ++i)
        r *= i;
    return r;
}

template<int d>
inline int multinomial(std::array<int, d> const& coeff) {
    int denom = 1, n = 0;
    for(int j=0; j<d; ++j) {
        n       += coeff[j];
        denom   *= factorial(coeff[j]);
    }
    return factorial(n) / denom;
}

//Index iterator for tensors
template<int R, int D>
struct SymmetricIndexIterator {

    static_assert(R >= 0, "Rank must be nonnegative");
    static_assert(D > 0, "Dimension must be greater than 0");

    //Constants
    enum {
        rank = R,
        dimension = D,
    };

    //Type macros
    typedef std::array<int, dimension>  DegreeIndex;
    typedef std::array<int, rank>       TensorIndex;
    
    //Constructors
    SymmetricIndexIterator() : index_(0) {
        for(int i=0; i<rank; ++i)
            seq_[i] = 0;
    }
    SymmetricIndexIterator(int pos) : index_(pos) { }
    SymmetricIndexIterator(DegreeIndex degrees) :
        seq_(degree_to_tensor(degrees)),
        index_(degree_to_array(degrees)) {}
    
    //Sequence comparison
    inline bool operator==(SymmetricIndexIterator const& other) const  { return index_ == other.index_; }
    inline bool operator!=(SymmetricIndexIterator const& other) const  { return index_ != other.index_; }
    inline bool operator<(SymmetricIndexIterator const& other) const   { return index_ <  other.index_; }
    inline bool operator<=(SymmetricIndexIterator const& other) const  { return index_ <= other.index_; }
    inline bool operator>(SymmetricIndexIterator const& other) const   { return index_ >  other.index_; }
    inline bool operator>=(SymmetricIndexIterator const& other) const  { return index_ >= other.index_; }
    
    //Iteration
    SymmetricIndexIterator& operator++() {
        for(int i=0; i<rank; ++i) {
            if(seq_[i] < dimension - 1) {
                ++seq_[i];
                ++index_;
                for(int j=i-1; j>=0; --j) {
                    seq_[j] = seq_[i];
                }
                return *this;
            }
        }
        index_ = -1;
        return *this;
    }
    
    //Postfix iteration
    inline SymmetricIndexIterator operator++(int) {
        auto tmp = *this;
        ++*this;
        return tmp;
    }
    
    //Accessors
    inline DegreeIndex        degree_index() const { return tensor_to_degree(seq_); }
    inline TensorIndex const& tensor_index() const { return seq_; }    
    inline int                array_index()  const { return index_; }
    
    //Helper methods
    static inline int degree_to_array(DegreeIndex const& degrees) {
        int index = 0, sum = 0;
        for(int i=dimension-1; i>0; --i) {
            sum += degrees[i];
            int w = 1;
            for(int j=0; j<dimension-i; ++j) {
                w *= sum + j;
            }
            index += w / factorial(dimension - i);
        }
        return index;
    }

    static inline TensorIndex degree_to_tensor(DegreeIndex const& degrees) {
        TensorIndex result;
        int j = 0;
        for(int i=0; i<dimension; ++i) {
            for(int k=0; k<degrees[i]; ++k) {
                result[j++] = i;
            }
        }
        return result;
    }

    static inline DegreeIndex tensor_to_degree(TensorIndex const& seq) {
        DegreeIndex result;
        for(int i=0; i<dimension; ++i) {
            result[i] = 0;
        }
        for(int i=0; i<rank; ++i) {
            ++result[seq[i]];
        }
        return result;
    }

private:
    TensorIndex seq_;
    int index_;
};

//The symmetric tensor
template<typename Scalar_t, int R, int D> 
struct SymmetricTensor {

    static_assert(R >= 0, "Rank must be nonnegative");
    static_assert(D > 0, "Dimension must be greater than 0");

    //Type macros
    typedef Scalar_t                                Scalar;
    typedef SymmetricIndexIterator<R, D> BaseIterator;
    typedef typename BaseIterator::DegreeIndex      DegreeIndex;
    typedef typename BaseIterator::TensorIndex      TensorIndex;
    
    //Constants
    enum {
        rank = BaseIterator::rank,
        dimension = BaseIterator::dimension,
        size = Binomial<rank + dimension - 1, rank>::value
    };
    
    //Iterator type
    friend class Iterator;
    struct Iterator : public BaseIterator {
        Iterator(BaseIterator const& base, SymmetricTensor* t) : 
            BaseIterator(base), owner(t) {  }
        inline Scalar& operator*() { return owner->coefficients[this->array_index()]; }
    private:
        SymmetricTensor* owner;
    };
    
    //Const iterator type
    friend class ConstIterator;
    struct ConstIterator : public BaseIterator {
        ConstIterator(BaseIterator const& base, const SymmetricTensor* t) : 
            BaseIterator(base), owner(t) {  }
        inline const Scalar& operator*() { return owner->coefficients[this->array_index()]; }
    private:
        const SymmetricTensor* owner;        
    };

    //Iterators
    inline Iterator begin()            { return Iterator(BaseIterator(), this); }
    inline Iterator end()              { return Iterator(BaseIterator(-1), this); }
    inline Iterator degree_index(DegreeIndex const& deg) { return Iterator(BaseIterator(deg), this); }
    inline Iterator tensor_index(TensorIndex const& i) { return Iterator(BaseIterator(BaseIterator::tensor_to_degree(i)), this); }
    inline ConstIterator begin() const { return ConstIterator(BaseIterator(), this); }
    inline ConstIterator end() const   { return ConstIterator(BaseIterator(-1), this); }
    inline ConstIterator degree_index(DegreeIndex const& deg) const { return ConstIterator(BaseIterator(deg), this); }
    inline ConstIterator tensor_index(TensorIndex const& i) const { return ConstIterator(BaseIterator(BaseIterator::tensor_to_degree(i)), this); }

    //Coordinate indexing (by degree)
    inline Scalar const& operator()(std::initializer_list<int> const& lst) const {
        DegreeIndex deg;
        std::copy(lst.begin(), lst.end(), deg.begin());
        return coefficients[BaseIterator::degree_to_array(deg)];
    }
    inline Scalar& operator()(std::initializer_list<int> const& lst) {
        DegreeIndex deg;
        std::copy(lst.begin(), lst.end(), deg.begin());
        return coefficients[BaseIterator::degree_to_array(deg)];
    }

    //Tensor addition
    SymmetricTensor operator+(SymmetricTensor const& other) const {
        SymmetricTensor result;
        for(int i=0; i<size; ++i) {
            result.coefficients[i] = 
                coefficients[i] + 
                other.coefficients[i];
        }
        return result;
    }
    
    //Scalar multiplication
    SymmetricTensor operator*(Scalar const& scalar) const {
        SymmetricTensor result;
        for(int i=0; i<size; ++i) {
            result.coefficients[i] = coefficients[i] * scalar;
        }
        return result;
    }
    
    //Tensor contraction (aka Einstein summation)
    template<int other_rank>
    SymmetricTensor<Scalar, rank - other_rank, dimension>
        inner_product(SymmetricTensor<Scalar, other_rank, dimension> const& other) const {
        
        enum { final_rank = rank - other_rank };
        static_assert(final_rank >= 0, "Invalid tensor size for index contraction.");
        typedef SymmetricTensor<Scalar, rank, dimension>        A_type;
        typedef SymmetricTensor<Scalar, other_rank, dimension>  B_type;
        typedef SymmetricTensor<Scalar, final_rank, dimension>  C_type;
        
        C_type result;
        for(auto iter=result.begin(); iter!=result.end(); ++iter) {
            Scalar w(0);
            for(auto iter2=other.begin(); iter2!=other.end(); ++iter2) {
                DegreeIndex a, b = iter2.degrees(), c = iter.degrees();
                for(int i=0; i<dimension; ++i) {
                    a[i] = c[i] + b[i];
                }
                w += Scalar(multinomial(b)) * (*iter2) * (*degree_index(a));
            }
            (*iter) = w;
        }
        return result;
    }
    
    //Symmetric tensor product (aka polynomial multiplication)
    template<int other_rank>
    SymmetricTensor<Scalar, other_rank + rank, dimension>
        outer_product(SymmetricTensor<Scalar, other_rank, dimension> const& other) const {
        
        enum { final_rank = rank + other_rank };
        typedef SymmetricTensor<Scalar, rank, dimension>        A_type;
        typedef SymmetricTensor<Scalar, other_rank, dimension>  B_type;
        typedef SymmetricTensor<Scalar, final_rank, dimension>  C_type;
        
        struct Visitor {
            A_type const& A;
            B_type const& B;
            DegreeIndex v, a, vpartial, apartial;
            Scalar w;
            
            Visitor(DegreeIndex const& vin, A_type const& Ain, B_type const& Bin) : 
                A(Ain), B(Bin), v(vin), w(0) {
                for(int i=0; i<dimension; ++i) {
                    a[i] = 0;
                    apartial[i] = A_type::rank;
                    if(i > 0) {
                        vpartial[dimension - i-1] = vpartial[dimension - (i+2)] + v[dimension - (i+1)];
                    }
                    else {
                        vpartial[dimension - i-1] = 0;
                    }
                }
            }
            
            void visit(const int k) {
            
                using namespace std;    //TODO: REMOVE ME
                
                
                if(k == dimension) {
                    DegreeIndex b;
                    cout << "Visiting: " << endl;
                    for(int i=0; i<dimension; ++i) {
                        b[i] = v[i] - a[i];
                        
                        cout << v[i] << '\t' << b[i] << '\t' << a[i] << endl;
                    }
                    auto a_val = Scalar(multinomial<dimension>(a)) * (*A.degree_index(a));
                    auto b_val = Scalar(multinomial<dimension>(b)) * (*B.degree_index(b));
                    
                    cout << "a_val = " << a_val << endl
                         << "b_val = " << b_val << endl;
                         
                    w += a_val * b_val;
                    return;
                }
                
                if(k > 0) {
                    apartial[k] = apartial[k-1] - a[k-1];
                } else {
                    apartial[k] = A_type::rank;
                }
                int l = max(0, apartial[k] - vpartial[k]), u = min(v[k], v[k] + apartial[k]);
                
                cout << "Looping k = " << k << ", lb = " << l << ", ub = " << u 
                     << ", apartial[k] = " << apartial[k] << ", vpartial[k] = " << vpartial[k] << endl;
                     
                     
                for(a[k] = l; a[k] <= u; ++a[k]) {
                    visit(k+1);
                }
            }
        };
        
        C_type result;
        for(auto iter=result.begin(); iter!=result.end(); ++iter) {
            auto deg = iter.degree_index();
            Visitor visitor(deg, *this, other);
            visitor.visit(0);
            *iter = visitor.w / Scalar(multinomial<dimension>(deg));
        }
        return result;
    }

private:
    //Coefficients
    Scalar coefficients[size];
};

//Prints out the components of a tensor as a polynomial in pretty latex formatting
template<typename Scalar, int rank, int dimension>
std::ostream& operator<<(std::ostream& os, SymmetricTensor<Scalar, rank, dimension> const& tensor) {
    bool need_plus = false;

    for(auto iter=tensor.begin(); iter!=tensor.end(); ++iter) {    

        auto deg = iter.degree_index();
        auto coef = Scalar(multinomial<dimension>(deg)) * (*iter);
        
        if(coef == 0) {
            continue;
        } else if(need_plus) {
            if(coef > 0) {
                os << " + ";
                if( coef != 1 ) os << coef;
            } else {
                os << " - ";
                if( coef != -1 ) os << -coef;
            }
        } else if(coef == -1) {
            os << '-';
        } else if(coef != 1) {
            os << coef;
        }
        need_plus = true;
        
        bool need_space = false;
        for(int i=0; i<dimension; ++i) {
            if(deg[i] == 0)
                continue;
            if(need_space) {
                os << ' ';
            }
            if(i >= 10) {
                os << "x_{" << i << '}';
            } else {
                os << "x_" << i;
            }
            if(deg[i] >= 10) {
                os << "^{" << deg[i] << '}';
            } else if(deg[i] > 1) {
                os << '^' << deg[i];
            }
            need_space = true;
        }
    }
    return os;
}

#endif

